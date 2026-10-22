# ============= Script for high-resolution N2 depth integration ==============
# This code depth-integrates N2 at the resolution of input bathymetry, using 
# an 'up-sampling' approach refined from Davies & Guinotte (2011).
# Here, a weighted shape-preserving spline is used on the local and all neighbouring
# cells to interpolate N2 at the depth of the seabed. Depth integration is carried
# out over the proportion of the water column specified by H_prop from the bottom
# up [H_prop = 0.1 - 10%,[H_prop = 1 - full water column], which may be capped to
# a maximum value [H_cap in metres] for deep waters.
# Options are depth integrating zN2 (multiplying N^2 by depth for body force), or 
# depth-averaging by (H*H_prop), i.e. for nearbed N2 calculation 
# output quality flags are produced to aid identifying areas of data paucity

required_packages <- c('ncdf4', 'terra', 'data.table')
suppressPackageStartupMessages({
  for (i in required_packages){
    library(i, character.only = TRUE)}
})

# =================== Input Settings ===================
# input 3D N2 file (coarse grid):
n2_nc_path   <- "02_secondary_data/NovaScotia/intermediary_files/summerN2.nc"

n2_varname   <- "N2"
# high-res bathymetry (here local LCC, meters; negative depths or positive down):
bathy_path   <- "clipped_bathy.tif"  
#local lcc proj4 definition:
lcc_crs <- paste("+proj=lcc +lat_0=46 +lon_0=-60 +lat_1=48 +lat_2=52",
                "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", sep = ' ')
# Outputs
out_dir <- r'(02_secondary_data\NovaScotia\intermediary_files)'
out_mean <- paste(out_dir,"N2_summer_100m_fullWC.tif", sep = '\\')
out_flag <-  paste(out_dir,"N2_winter_100m_fullWC.tif", sep = '\\') 
# categorical flag containing 4 classes:
# 0 - local N2 water column includes seafloor
# 1 - ngb N2 spline includes seafloor, max depth values based on ngb N2
# 2 - reduced window (==H_prop) due to insufficient data in ngb at max depth
# 3 - no data within H_prop, returned NA
# These flags may help tune H_prop below to a sensible value
out_qflag <-  paste(out_dir,"N2_summer_100m_fullWC.tif", sep = '\\')
# quantitative quality flag
# percentage of theoretical values used in nearbed range across 9 cell ngb

# integration window, change to desired nearbed fraction, H_prop=1 <- whole column
H_prop <- 1
H_cap  <- 6000 # effectively removes cap as water depth <6000 everywhere

# neighborhood‑weighted fitting
w_local      <- 1.0 #local cell weight for spline fits
w_nei        <- 0.25 #8 neighbour weight for smoothing_splines
# hermitic splines are impossible due to non-monotone nature of input data,
# instead control shape of splines by constraining smoothing
spar_in <- list(low=0, high=0.2) 
# optional (and experimental, keep = FALSE unless you know what you are doing)
use_winsor   <- FALSE
winsor_probs <- c(0.005, 0.995)

min_pts_main  <- 6
min_pts_local <- 6

# optional smoothing, this os experimental and not used:
# do_smooth       <- FALSE
# smooth_sigma_px <- 0.8
# smooth_radius   <- 2L

# WOA mid‑depth levels (positive down via abs())
mid_levels <- depth_levels[-length(depth_levels)]+diff(depth_levels)
z_mid <- abs(mid_levels)

#intergration options
depth_factor <- FALSE    # multiplies N2 by respective N2 depth within intergrand
depth_avg <- TRUE        # divides by (H*H_prop) for depth-averaged integration  

# =================== HELPERS ===================

winsorize <- function(y, probs=c(0.005,0.995)) {
  qs <- stats::quantile(y, probs=probs, na.rm=TRUE)
  pmax(qs[1], pmin(qs[2], y))
}

fit_spline <- function(z, y, w=NULL) {
  o <- order(z)
  if (is.null(w)) {
    stats::smooth.spline(x=z[o], y=y[o], control.spar=spar_in)
  } else {
    stats::smooth.spline(x=z[o], y=y[o], w=w[o], control.spar=spar_in)
  }
}

# 8‑point Gauss–Legendre integral (with optional interval average)
GL_rules <- pracma::gaussLegendre(8,-1,1) #8 should do, but can be changed up
.gl8x <- GL_rules$x
.gl8w <- GL_rules$w

int_over_interval <- function(f, a, b) {
  if (!is.finite(a) || !is.finite(b) || b <= a) return(NA_real_)
  zz <- ((b - a)/2) * .gl8x + (a + b)/2
  yy <- f(zz)
  if (all(!is.finite(yy))) return(NA_real_)
  #change of interval to [a,b], generalised Gauss-Legendre integral:
  sum(.gl8w * yy) * (b - a)/2 / 
    ifelse(depth_avg,(b - a),1) #depth_avg division if required
}

# 8‑neighbor IDs
ngb8_ids <- function(ix, iy, nx, ny) {
  ii <- pmax(1, pmin(nx, ix + c(-1,0,1)))
  jj <- pmax(1, pmin(ny, iy + c(-1,0,1)))
  g  <- as.data.table(expand.grid(ix = ii, iy = jj))
  g  <- g[!(ix == !!ix & iy == !!iy)]
  g[, id := (iy-1)*nx + ix]
  unique(g$id)
}
# for [optional & experimental] smoothing
gauss_kernel <- function(sigma_cells=0.8, radius=2L) {
  r <- max(1L, as.integer(radius))
  x <- -r:r; y <- -r:r
  K <- outer(x, y, function(a,b) exp(-(a*a + b*b)/(2*sigma_cells^2)))
  K / sum(K)
}

# =================== READ N2 ===================

nc <- nc_open(n2_nc_path)
on.exit(try(nc_close(nc), silent=TRUE), add=TRUE)

var <- nc$var[[n2_varname]]
if (is.null(var)) stop("Variable '", n2_varname, "' not found in ", n2_nc_path)

dimnames <- sapply(var$dim, `[[`, "name")
dimsz    <- sapply(var$dim, `[[`, "len")

z_dim_idx <- which.max(dimsz == length(z_mid))
if (length(z_dim_idx)!=1) stop("Cannot identify depth dimension")

xy_idx <- setdiff(seq_along(dimsz), z_dim_idx)
nx <- dimsz[xy_idx[1]]; ny <- dimsz[xy_idx[2]]; nz <- dimsz[z_dim_idx]
N2_raw <- ncvar_get(nc, n2_varname)

# permute to [ny, nx, nz]
perm <- c(NA,NA,NA); perm[z_dim_idx] <- 3; perm[xy_idx] <- c(2,1)
N2_yxz <- aperm(N2_raw, perm=perm)

get1d <- function(names) {
  for (nm in names) {
    v <- try(ncvar_get(nc, nm), silent=TRUE)
    if (!inherits(v,"try-error") && !is.null(v)) return(as.numeric(v))
  }
  NULL
}
lon <- get1d(c("lon","longitude","x"))
lat <- get1d(c("lat","latitude","y"))
if (length(lon)!=nx || length(lat)!=ny) stop("Unexpected lon/lat dims")

# ensure row 1 is the northmost row
if (lat[1] < lat[ny]) { N2_yxz <- N2_yxz[ny:1,,,drop=FALSE]; lat <- rev(lat) }

# =================== NEIGHBORS ===================

neighbors_list <- vector("list", ny*nx)
for (iy in 1:ny) for (ix in 1:nx) {
  neighbors_list[[(iy-1)*nx + ix]] <- ngb8_ids(ix, iy, nx, ny)
}

# ================= FIT SPLINES & Count Valid input for qflags ==================

cat("Fitting neighborhood-weighted splines and counting valid samples per depth ...\n")

cell_model <- vector("list", ny*nx)

for (iy in 1:ny) for (ix in 1:nx) {
  id <- (iy-1)*nx + ix
  
  # collect 3×3 column values (central + 8 neighbors)
  ngb <- c(id, neighbors_list[[id]])
  
  # --- counts9 over depth (max 9 per z level) ---  
  counts9 <- integer(nz)
  for (nid in ngb) {
    jy <- ((nid-1) %/% nx) + 1
    jx <- ((nid-1) %%  nx) + 1
    yv <- N2_yxz[jy, jx, ]
    counts9 <- counts9 + as.integer(is.finite(yv))
  }
  # clamp to [0,9]
  counts9[counts9 < 0] <- 0L; counts9[counts9 > 9] <- 9L
  
  # build MAIN spline (weighted neighborhood)
  z_all <- z_mid
  z_vec <- c(); y_vec <- c(); w_vec <- c()
  
  # central first
  yC <- N2_yxz[iy, ix, ]
  okC <- is.finite(yC)
  if (any(okC)) {
    zc <- z_all[okC]; yc <- yC[okC]
    if (use_winsor) yc <- winsorize(yc, winsor_probs)
    z_vec <- c(z_vec, zc); y_vec <- c(y_vec, yc); w_vec <- c(w_vec, rep(w_local, length(zc)))
  }
  
  # neighbors
  for (nid in neighbors_list[[id]]) {
    jy <- ((nid-1) %/% nx) + 1
    jx <- ((nid-1) %%  nx) + 1
    yN <- N2_yxz[jy, jx, ]
    okN <- is.finite(yN)
    if (any(okN)) {
      zn <- z_all[okN]; yn <- yN[okN]
      if (use_winsor) yn <- winsorize(yn, winsor_probs)
      z_vec <- c(z_vec, zn); y_vec <- c(y_vec, yn); w_vec <- c(w_vec, rep(w_nei, length(zn)))
    }
  }
  
  fun_main <- NULL; zmin_main <- NA_real_; zmax_main <- NA_real_
  if (sum(is.finite(y_vec)) >= min_pts_main) {
    # multiplication of depth by N2 for body force when depth_factor = TRUE
    y_vec <- ifelse(depth_factor,z_vec,1) * y_vec 
    fitM <- try(fit_spline(z_vec, y_vec, w_vec), silent=TRUE)
    if (!inherits(fitM,"try-error")) {
      fun_main  <- function(zz) as.numeric(stats::predict(fitM, zz)$y)
      zmin_main <- min(z_vec); zmax_main <- max(z_vec)
    }
  }
  
  # local-only spline (for process flag=0)
  fun_local <- NULL; zmin_local <- NA_real_; zmax_local <- NA_real_
  if (sum(is.finite(yC)) >= min_pts_local) {
    fitL <- try(fit_spline(z_all[okC], yC[okC]), silent=TRUE)
    if (!inherits(fitL,"try-error")) {
      fun_local  <- function(zz) as.numeric(stats::predict(fitL, zz)$y)
      zmin_local <- min(z_all[okC]); zmax_local <- max(z_all[okC])
    }
  }
  
  cell_model[[id]] <- list(
    fun_main   = fun_main,  zmin_main = zmin_main, zmax_main = zmax_main,
    fun_local  = fun_local, zmin_local= zmin_local, zmax_local= zmax_local,
    counts9    = counts9    # per-depth valid count across 3×3
  )
}

# =================== Load BATHY ===================================

terraOptions(memfrac=0.3)
bathy <- -rast(bathy_path)
if (is.na(crs(bathy))) crs(bathy) <- lcc_crs
if (crs(bathy) != lcc_crs) bathy <- project(bathy, lcc_crs, method="bilinear")
# only run the following if it's not arleady in memory, 
# they are the same across all runs
if(!exists('r_id_lcc') | 
   inherits(try(str(r_id_lcc@pntr), silent = T), 'try-error') || 
   crs(r_id_lcc)!=crs(bathy)){
  # QA rasters: 0-based row/col/id of WOA grid, projected tobathy (nearest)
  r_ll <- rast(ncols=nx, nrows=ny,
               extent=ext(min(lon), max(lon), min(lat), max(lat)),
               crs="EPSG:4326")
  m_id  <- matrix(0:(ny*nx-1),             nrow=ny, ncol=nx, byrow=TRUE)
  # m_row <- matrix(rep(0:(ny-1), each=nx),  nrow=ny, ncol=nx, byrow=TRUE)
  # m_col <- matrix(rep(0:(nx-1), times=ny), nrow=ny, ncol=nx, byrow=TRUE)
  # 
  r_id_ll   <- r_ll; values(r_id_ll)   <- m_id
  # r_row0_ll <- r_ll; values(r_row0_ll) <- m_row
  # r_col0_ll <- r_ll; values(r_col0_ll) <- m_col
  # 
  r_id_lcc   <- project(r_id_ll,   bathy, method="near")
  # r_row0_lcc <- project(r_row0_ll, bathy, method="near")
  # r_col0_lcc <- project(r_col0_ll, bathy, method="near")
  # 
  # writeRaster(r_row0_lcc, "n2_row0.tif", overwrite=TRUE)
  # writeRaster(r_col0_lcc, "n2_col0.tif", overwrite=TRUE)
  writeRaster(r_id_lcc,   "n2_id0.tif",  overwrite=TRUE)
}
# =================== STREAMED COMPUTE & WRITE ===================

out_val <- rast(bathy); names(out_val) <- "N2_nearbed"
out_flg <- rast(bathy); names(out_flg) <- "flags"
out_qf  <- rast(bathy); names(out_qf)  <- "qflag"    

b1 <- writeStart(out_val, filename=out_mean, overwrite=TRUE,
                 filetype="GTiff", gdal=c("COMPRESS=DEFLATE","ZLEVEL=6","TILED=YES"))
b2 <- writeStart(out_flg, filename=out_flag, overwrite=TRUE,
                 filetype="GTiff", gdal=c("COMPRESS=DEFLATE","ZLEVEL=6","TILED=YES"))
b3 <- writeStart(out_qf,  filename=out_qflag, overwrite=TRUE,
                 filetype="GTiff", gdal=c("COMPRESS=DEFLATE","ZLEVEL=6","TILED=YES"))

b <- b1
cat("Flags: 0= full with local throughout profile, 1=neighbor-assisted full, 2=reduced window, 3=NA\n")
cat("qflag: % of valid (3×3) WOA samples available in [z-H, z]\n")

for (i in 1:b$n) {
  r0 <- b$row[i]; nr <- b$nrows[i]
  zb  <- readValues(bathy,    r0, nr)   # depth
  pid <- readValues(r_id_lcc, r0, nr)   # 0-based parent id
  
  if (all(zb <= 0, na.rm=TRUE)) zb <- abs(zb)
  
  v_mean <- rep(NA_real_, length(zb))
  v_flag <- rep(3L,        length(zb))
  v_qflg <- rep(NA_real_,  length(zb))   # 
  
  for (k in seq_along(zb)) {
    z   <- zb[k]; id0 <- pid[k]
    if (!is.finite(z) || is.na(id0)) { v_flag[k] <- 3L; next }
    cid <- id0 + 1L
    cm  <- cell_model[[cid]]
    if (is.null(cm) || is.null(cm$fun_main)) { v_flag[k] <- 3L; next }
    
    # near-bed window, or full water column if H_prop=1 and H_cap>max_depth
    H     <- min(H_prop * z, H_cap)
    # trim top of window to z_mid[1] for full WC (rather than 0), to avoid flag
    # defaulting to reduced window
    a_req <- ifelse(z - H < z_mid[1],z_mid[1],z - H)
    b_req <- z
    
    # qflag computation (based on raw availability of ngb cells, not weights)
    # indices of WOA mid-levels inside [a_req, b_req]
    idx_win <- which(z_mid >= a_req & z_mid <= b_req)
    n_theory <- 9L * length(idx_win)
    if (n_theory > 0L) {
      n_valid <- sum(cm$counts9[idx_win], na.rm=TRUE)
      v_qflg[k] <- 100.0 * (n_valid / n_theory)
    } else {
      v_qflg[k] <- NA_real_
    }
    #BUG FIX needed for H_prop=1, a_req will be 0, defaulting to reduced window
    
    # process flag + mean
    if (!is.null(cm$fun_local) &&
        a_req >= cm$zmin_local && b_req <= cm$zmax_local) {
      # here, fun_local could theoretically be used
      v_mean[k] <- int_over_interval(cm$fun_main, a_req, b_req)
      v_flag[k] <- 0L
      next
    }
    
    if (a_req >= cm$zmin_main && b_req <= cm$zmax_main) {
      v_mean[k] <- int_over_interval(cm$fun_main, a_req, b_req)
      v_flag[k] <- 1L
      next
    }
    
    # reduced window with main
    a_eff <- max(a_req, cm$zmin_main)
    b_eff <- min(b_req, cm$zmax_main)
    if (b_eff > a_eff) {
      v_mean[k] <- int_over_interval(cm$fun_main, a_eff, b_eff)
      v_flag[k] <- 2L
    } else {
      v_flag[k] <- 3L
    }
  }
  
  writeValues(out_val, v_mean, r0, nr)
  writeValues(out_flg, v_flag, r0, nr)
  writeValues(out_qf,  v_qflg, r0, nr)
}

writeStop(out_val)
writeStop(out_flg)
writeStop(out_qf)

cat("Finished main outputs:\n  - ", out_mean, "\n  - ", out_flag, "\n  - ", out_qflag, "\n", sep="")

# optional smoothing of the mean field
# if (isTRUE(do_smooth)) {
#   cat("Writing smoothed copy ...\n")
#   K <- gauss_kernel(sigma_cells=smooth_sigma_px, radius=smooth_radius)
#   fn_out <- paste(strsplit(out_mean,'\\.')[[1]][1],"smooth.tif",sep = "_")
#   N2_smooth <- focal(rast(out_mean), w=K, fun=sum, na.policy="omit",
#                      filename=fn_out, overwrite=TRUE)
#   cat("Finished optional outputs:\n",fn_out)
# }

cat("QA rasters:\n  - n2_row0.tif\n  - n2_col0.tif\n  - n2_id0.tif\n")
cat("Done.\n")
