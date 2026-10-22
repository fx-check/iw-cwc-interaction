
## Script to integrate uv data with DEM tiles
## Rationale: expand time dimension to 365 for each DEM tile
#  Here we write new .nc files, one with t=365 and one with mean annual E flux,  
# from numeric integration of spectral power density across wave number bins, 
# and daily current data as per Nikurashin & Ferrari (2011)
# We avoid writing each daily flow aligned wave number matrix to save disk space, 
# the manipulation of the raw matrix is fast in memory (on the fly)
### Method outline
## Grab u and v from 2022_uv.nc file and commit appropriate cell data to environment
## with the DEM tile names/data. Write vector of u and v to list variables
## Use this data to flow-align wave number matrix and mask spectral power matrix
## integrate spectral power density and by summing across binned wave numbers
## write the summed bins to environment. The workflow is vectorised and parallelised

required_packages <- c('ncdf4','stringr','lubridate','future.apply','parallelly',
                       'stars')
for (i in required_packages){
suppressPackageStartupMessages(library(i, character.only = T))
}
t_var <- as.character(as_date('2021-12-31')+days(1:365))

#
# prior to heavy computations, it may be worth removing large vars from memory
tail(sort(sapply(.GlobalEnv, object.size)), 10)
# rm(list = names(tail(sort(sapply(.GlobalEnv, object.size)), 10)[1:8]))
# Preliminaries -----------------------------------------------------------

stopifnot("Expected psd dimensions are 128 by 128" = (unique_psd_dims==128))
fft_path <- r'(C:\Users\fbutschek\Desktop\hycom\fft_output)'
  # '03_derived_data/Blake/fft_tiles/'
# Quick loop to double check our spectral dims, this allows precomputing ky & kx:
for(i in ls(fft_env)){
  cur_var <- get(i, envir = fft_env)
  if(!all.equal(cur_var$psd_dim,c(128,128))){
    print(paste(i, 'has dimensions', cur_var$psd_dim[1], ',', cur_var$psd_dim[1]))}
}
# all spectral wave vectors are the same, so we can precompute kx and ky
# for the unique dimensions


# Precompute wave vectors -------------------------------------------------
L <- unique_psd_dims[1]; dr <- 35
dfr <- 1/ (L*dr) # frequency in cycles per unit (dr)
dl = dk = 2 * pi *dfr # wave vector spacing, identical along and across flow

# replicate in wave vectors below
k_xy <- (-L/2):(L/2 - 1) * (dk) # multiplication of index by angular wavenumber spacing
#matrices from k_x and k_y, :
KX =  matrix(rep(k_xy, each = L), nrow = L, ncol = L)
KY = t(KX)
# sanity check: 
# identical(f_out$ky,KY); identical(f_out$kx,KX) #check passes against last f_out

# set up variables with path names to grab nc files from within parallel loop -----------------


uv_path <- r'(C:\Users\fbutschek\Desktop\hycom\2022\u_v_bottom\uv_geos_utm_near_revised.nc)'
N2_path=r'(C:\Users\fbutschek\Desktop\hycom\2022\seasonal_bottom_N2_data_UTM.nc)'
gamma_path <- r'(C:\Users\fbutschek\Desktop\hycom\2022\gamma_grid.nc)'
uv_geos_2022_UTM18N <- read_stars(uv_path)
# Add grid convergence offset to fft_env ----------------------------------
# if you haven't added this previously
gamma_nc <- nc_open(gamma_path)

gamma_grid <- ncvar_get(gamma_nc, 'gamma')

nc_close(gamma_nc)
for (i in ls(fft_env)){
  lst <- get(i, fft_env)
  lst$gamma<- gamma_grid[lst$nx, lst$ny]
  assign(i, lst, fft_env)
}
rm(gamma_grid)
# recover latitudes for Coriolis frequency --------------------------------

# the nc files themselves need to be opened from within the worker to allow 
# parallel sessions accessing the data
northings <- st_get_dimension_values(uv_geos_2022_UTM18N, 
                                     which(names(
                                       st_dimensions(uv_geos_2022_UTM18N))=='y'))
eastings <- st_get_dimension_values(uv_geos_2022_UTM18N, 
                                     which(names(
                                       st_dimensions(uv_geos_2022_UTM18N))=='x'))
#t=1 of N2 raster keeps the follwoing lighter:
coord_star <- st_dim_to_attr(read_stars(N2_path), which(names(
  st_dimensions(uv_geos_2022_UTM18N))== 'x'))
coord_star$y <- st_dim_to_attr(coord_star, which(names(
  st_dimensions(uv_geos_2022_UTM18N))== 'y'))

geo_coords <- oce::utm2lonlat(easting = coord_star$x,
                              northing = coord_star$y, zone = 18)
lats <- geo_coords$latitude[22,] #for coriolis only, mid long sufficient
rm(uv_geos_2022_UTM18N); rm(geo_coords)
# longs <- geo_coords$longitude[,1]

# Parallelised worfklow ---------------------------------------------------
kxy_2D <- cbind(as.vector(KX), as.vector(KY))              # (L*L x 2)

#for N2, we'll expand the 2-d data to a time dim using the following:
spring_end <- 79 
summer_end <- 263
# min and max wave vectors exceed viable lee wave bands, so let's stick with
# bins based on non-rotated space and add +/- infinity 
across_flow_bin_edges <- seq(k_xy[1], k_xy[length(k_xy)], length.out = L)
k_mid_bins <- head(across_flow_bin_edges,-1)+ abs(head(across_flow_bin_edges,-1)-
                                               tail(across_flow_bin_edges,-1))/2
across_flow_bin_edges <- c(-Inf, across_flow_bin_edges, Inf) #avoids NAs 
across_flow_bins <- cut(seq(-L/2*dk, L/2*dk, length.out = L+1), 
                        across_flow_bin_edges) #for matching

#extend length to match, exact position inmaterial as will always be masked:
k_mid_bins <- c(-2*across_flow_bin_edges[2], #beyond negative extent (avoids -Inf)
                k_mid_bins, #mid bins from differencing, extend by 1 assuming same delta
                across_flow_bin_edges[L]+(across_flow_bin_edges[L]-across_flow_bin_edges[L-1])/2,
                2*across_flow_bin_edges[L]) #beyond positive extent (avoids +Inf)

k <- length(across_flow_bin_edges) # as -INF to INF included, no subtraction needed
stopifnot(length(k_mid_bins)==k)


# nc open separate copies for workers ---------------------------------------------------
# this function prevents file handling issues by creating a local copy for each 
# worker, which is desireable if the .nc files are moderate in size and the
# per-worker time low, preventing simultaneous opening of nc files on disk by 
# different workers. 
.local_nc_cache <- new.env(parent = emptyenv())

with_local_nc <- function(path) {
  # Use original absolute path as key (unique even if basenames collide)
  key <- normalizePath(path, winslash = "/", mustWork = TRUE)
  
  # If not cached for this worker, copy it
  if (is.null(.local_nc_cache[[key]])) {
    pid <- Sys.getpid()
    tmpd <- file.path(tempdir(), sprintf("nc_worker_%d", pid))
    dir.create(tmpd, showWarnings = FALSE, recursive = TRUE)
    
    tpath <- file.path(tmpd, basename(path))
    # Avoid overwriting if it already exists
    if (!file.exists(tpath)) file.copy(path, tpath, overwrite = FALSE)
    
    .local_nc_cache[[key]] <- tpath
  }
  
  .local_nc_cache[[key]]
}

## Worker function: process one tile ---------------------------------------

process_tile <- function(curr_fft,eapp = T, local_copies = T,...) {
  if (eapp ==FALSE) {#for tests or other apply fx,
    curr_fft <- get(curr_fft, envir = fft_env)}
  ### load all u and v timesteps for desired cell per worker-------------------
  uv_path_eff <- if (isTRUE(local_copies)) with_local_nc(uv_path) else uv_path
  uv_geos_utm <- nc_open(uv_path_eff)
  on.exit(try(nc_close(uv_geos_utm),silent=T), add = T)
  vars <- names(uv_geos_utm$var)
  if (!all(c("geostrophic_u_bottom","geostrophic_v_bottom") %in% vars)) {
    stop(sprintf("Expected velocity variables not found. Opened file was: %s",
                 uv_geos_utm$filename))
  }
  # v_geos_utm <- nc_open(filename = v_path)
  N2_path_eff <- if (isTRUE(local_copies)) with_local_nc(N2_path) else N2_path
  seasonal_N2 <- nc_open(N2_path_eff)
  on.exit(try(nc_close(seasonal_N2),silent=T), add = T)
  vars <- names(seasonal_N2$var)
  if (!all(c("winter_N2","summer_N2") %in% vars)) {
    stop(sprintf("Expected N2 variables not found. Opened file was: %s",
                 seasonal_N2$filename))
  }
  u_geos <- ncvar_get(uv_geos_utm, 'geostrophic_u_bottom', 
                               start = c(curr_fft$nx, curr_fft$ny, 1), 
                               count = c(1,1,-1)) 
  v_geos <- ncvar_get(uv_geos_utm, 'geostrophic_v_bottom', 
                               start = c(curr_fft$nx, curr_fft$ny, 1), 
                               count = c(1,1,-1)) 
  # nc_close(uv_geos_utm)
  theta <- atan2(v_geos,u_geos) #flow direction from east, CCW
  gamma <- curr_fft$gamma #grid convergence angle
  theta_true <- theta + gamma
  cos_theta <- cos(theta_true)
  sin_theta <- sin(theta_true)
  
  ### load N2 data for each timestep for the desired cell ---------------------
  
  N2_win <- ncvar_get(seasonal_N2, 'winter_N2', 
                                   start = c(curr_fft$nx, curr_fft$ny), 
                                   count = c(1,1))
  N2_sum <- ncvar_get(seasonal_N2, 'summer_N2', 
                                   start = c(curr_fft$nx, curr_fft$ny), 
                                   count = c(1,1))
  # nc_close(seasonal_N2) #interferese with on.exit() despite try()
  N2 <- c(rep(N2_win,spring_end),
                   rep(N2_sum,summer_end-spring_end),
                   rep(N2_win,365-summer_end)
  ) # vector with length 365
  ### --- compute frequency thresholds for the desired cell ---------------------
  u_mag <- sqrt(u_geos^2 +v_geos^2)
  k_max <- sqrt(N2)/u_mag # = N/U
  # coriolis and low frequency thresholds
 
  lat <- lats[curr_fft$ny]
  f_w<- oce::coriolis(lat)
  
  k_min <- abs(f_w)/u_mag # |f|/|U|, abs(f) is important on SH
  #safeguard against divisions by 0: 
  valid_flow <- is.finite(u_mag) & u_mag > 0 &
    is.finite(N2) & N2 > 0
  
  ### --- rotate k to along-flow for every t (matrix multiply) ----
  # set up 3d efficient vectorised workflow:
  coeff_al <- rbind(cos_theta,  sin_theta)   # coefficient for k_along
  k_along_2d <- kxy_2D %*% coeff_al #along flow aligned wavenumbers, with rows 
  #representing vectorised cells and columns equal to time dimension
  ### load and rearrange fft results for  integration --------------------------
  ## load tile's power spectrum from fft
  nc_fft <- nc_open(paste(fft_path,curr_fft$filename, sep = '/'))
  on.exit(try(nc_close(nc_fft),silent = T), add = T)
  #to convert spectral power matrix to spectral density, multiply by area
  psd <- ncvar_get(nc_fft, 'spectral_power_matrix') *(L*dr)^2 
  psd_vec <- as.vector(psd) #vectorise psd as well as kmag for integral
  # k_al <- rep(k_mid_bins[-1,-length(k_mid_bins)],L) 
  kmag_vec <-  as.vector(sqrt(KX^2+KY^2)) #vectorised cells, same as k_along_2d
  # Frequency windows per day,  |f|/|U| < k < |N|/|U|  ---------------
  min_bins <- cut(k_min, across_flow_bin_edges)
  kmin_idx <- match(min_bins,across_flow_bins) 
  max_bins <- cut(k_max, across_flow_bin_edges)
  kmax_idx <- match(max_bins,across_flow_bins) 
  #indices to be EXCLUDED, thus +/- 1 or min/max respectively, +2 to account -kmax_idx
  boundries <- data.frame('kmin_idx_neg' = k - kmax_idx + 2, 
                          'kmax_idx_neg'= (k%/%2+1)-(kmin_idx %% (k%/%2))-1,
                          'kmin_idx_pos' = kmin_idx+1, 'kmax_idx_pos'= kmax_idx-1)
  keep_idx <- apply(boundries, MARGIN = 1, FUN = function(x){
    if(all(sort(as.vector(unlist(x)))==as.vector(unlist(x)))){
      c(seq.int(x[1],x[2]),seq.int(x[3],x[4]))}})#mask of indices to keep

# along-flow integration Froude numbers & energy flux ------------------------------------
  #for all lists (except bins) in the across_flow_integral, sum using keep_idx
  # pre-allocate vectors
  E_flux = E_flux_blocked = numeric(length(t_var))
  Fr_l = linear_wave_penalty = rep(NA_real_, length(t_var))
  # |k|/\mathbf{|k|}P(k,l)dl across flow integration to be implemented first:
  spectral_weight <- abs(k_along_2d) / kmag_vec
  spectral_weight[kmag_vec==0, ] <- 0 # replace NA from 0 division above
  prod_mat <- spectral_weight * psd_vec   # (L^2 x 365), vectorized
  
  for (j in seq_len(length(t_var))) {
    idx <- findInterval(k_along_2d[, j], across_flow_bin_edges)   # bin index per point
    # mask to propagating band:
    keep <- idx %in% keep_idx[[j]]
    if (!any(keep)) next
    # integral over P(k,l) for characteristic topographic height and Froude no ----
    char_topo <- sqrt( #convert PSD to P_DFT for mean square amplitude: 
      sum(psd_vec[keep]/(L*dr)^2, na.rm = T) # discrete integral as in Perron08
    )
    Fr_l[j] <- sqrt(N2[j])*char_topo/u_mag[j]; Fr_c <- 0.4 #lee wave and critical Fr numbers
    linear_wave_penalty[j] <- ifelse( Fr_l[j] <= Fr_c,1, (Fr_c/ Fr_l[j])^2) 
    # sum by bin without reshape, this returns a 1d matrix per timestep with ixky in dimnames:
    acc <- rowsum(prod_mat[keep, j, drop = FALSE], group = idx[keep]) * dl
    eff_top <- rep(0, k) #zero-pads to k, with 0 for out-of bound k (beyond f/U or N/U)
    eff_top[as.integer(rownames(acc))] <- acc[, 1] / (2*pi)
    U2 <- u_mag[j]^2;  f2 <- f_w^2
    k2 <- k_mid_bins[keep_idx[[j]]]^2
    along_flow_integrand <- eff_top[keep_idx[[j]]] * 
      # sqrt(N2[j] - U2*k2) * sqrt(U2*k2 - f2) #for testing
      sqrt(pmax(0, N2[j] - U2*k2)) * sqrt(pmax(0, U2*k2 - f2)) #safeguard
    if (!valid_flow[j]) {
      E_flux[j] <- NA_real_
      E_flux_blocked[j] <- NA_real_
      next
    }
    E_flux[j] <- 1025 * u_mag[j] / (2*pi) * sum(along_flow_integrand)*dk
    E_flux_blocked[j] <- E_flux[j] * linear_wave_penalty[j]
  }
  
  
  lst <- list('u_geos' = u_geos, 'v_geos' = v_geos, 'speed' = u_mag, 
              'k_min' = k_min, 'k_max' = k_max,'coriolis_frequency' = f_w, 
              'winter_N2' = N2_win, 'summer_N2' = N2_sum, 'Fr_L' = Fr_l,
              'lin_wave_penalty' = linear_wave_penalty)
  lst$E_flux <- E_flux*1000 #in mW/m2
  lst$E_flux_blocked <- E_flux_blocked*1000 # accounts for upstream blocking
  # return(assign(i, curr_fft, envir = fft_env))
  
  lst
}

# the function was tested extensively, it seems to work fine
# peak RAM use is about 400-500MB, on this machine with 32GB RAM, 22 parallel cores
# should not max out RAM, but fewer might make it more stable


# Implement function using parallel processing ----------------------------


# explicitly export the globals that the worker needs, this can lead to hick-ups
# thus commented out in the function
# globals_to_send <- c("fft_path", "uv_path", "N2_path", "t_var",
#                      "KX", "KY", "kxy_2D", "kmag_vec",
#                      "across_flow_bin_edges","across_flow_bins", "k_mid_bins", 
#                      "dk", "k", "spring_end", "summer_end")



# Remove all-na cells from fft --------------------------------------------
# some cells with u,v N2 NA values are introduced as a result of reprojection
# the associated fft tile are eliminated here
uv_geos_utm <- nc_open(filename = uv_path)
NA_cells <- c()
for(i in ls(fft_env)){
  curr_fft <- get(i, envir = fft_env) #for test or other apply fx, comment out for eapply
  ### load all u and v timesteps for desired cell per worker-------------------
  u_geos <- ncvar_get(uv_geos_utm, 'geostrophic_u_bottom', 
                      start = c(curr_fft$nx, curr_fft$ny, 1), 
                      count = c(1,1,-1)
                      )
  if(all(is.na(u_geos))){
    rm(list = i,envir = fft_env)
    NA_cells <- c(NA_cells,i)
  }
}
nc_close(uv_geos_utm)
length(NA_cells)
## TESTS -------------------------------------------------------------------
j <- sample.int(length(ls(fft_env)),1)
i <- ls(fft_env)[j]
curr_fft <- get(i, envir = fft_env)
# process_tile(curr_fft) #uncomment to check output

peakRAM::peakRAM(process_tile(curr_fft, local_copies = F))
options(parallelly.debug = TRUE, future.debug = TRUE)
options(parallelly.makeNodePSOCK.setup_strategy = "sequential")
options(parallelly.makeNodePSOCK.validate = FALSE)
plan(multisession, workers = 4)
test_env <- new.env()
for (i in 1:4){
  nm <- ls(fft_env)[i]
  assign(nm, get(nm, fft_env), test_env)}
out <- future_eapply(test_env, FUN = function(i) process_tile(i),
                     # future.packages = c("ncdf4","data.table","lubridate","oce"),
                     future.seed = TRUE)

# out <- future_vapply(ls(test_env), FUN.VALUE = lst,
#                      FUN = function(i) process_tile(i, eapp = F),
#                      # future.packages = c("ncdf4","data.table","lubridate","oce"),
#                      future.seed = TRUE)
## ROLL-OUT ----------------------------------------------------------------
options(future.debug = FALSE, parallelly.debug = FALSE)
options(parallelly.makeNodePSOCK.setup_strategy = "sequential")
options(parallelly.makeNodePSOCK.validate = FALSE)
# roll out across all tiles:
workers <- min(8,max(1L, parallel::detectCores() - 2L))
plan(multisession, workers = workers)

options(future.debug = F, parallelly.debug = F)

out_ls3 <- future_eapply(
  env = fft_env,
  FUN = function(i) process_tile(i),
  # future.globals = c(globals_to_send, "process_tile"),
  # future.packages = c('ncdf4','data.table', 'lubridate','oce'),
  future.seed = TRUE)

# out_ls <- future_vapply(
#   ls(fft_env), FUN.VALUE = lst,
#   FUN = function(i) process_tile(i),
#   # future.globals = c(globals_to_send, "process_tile"),
#   future.packages = c('ncdf4','data.table', 'lubridate','oce'),
#   future.seed = TRUE)

#results:
all_E_mean <- vapply(out_ls3, 
                     FUN = function(x){log10(mean(x$E_flux, na.rm = T))}, 
                     numeric(1))
all_E_mean_blocked <- vapply(out_ls3, 
                             FUN = function(x){log10(
                               mean(x$E_flux_blocked, na.rm = T))}, 
                             numeric(1))
#ratios:
range(10^all_E_mean_blocked/10^all_E_mean); mean(10^all_E_mean_blocked/10^all_E_mean)
hist(all_E_mean_blocked)

# Write mean energy flux to nc file for mapping ---------------------------


if(!exists("date_dim")){
  date_dim <- ncdim_def(name = 'date', units = 'days since 2022-01-01', 
                      vals = 0:364, unlim = F)}
# double resolution for x dimension to account for tiles, this is a bit hacky but it works
x_vals <- eastings
hr_xval <- seq(x_vals[1], by = unique(diff(x_vals))/2,length.out = 2*length(x_vals))
all(x_vals==hr_xval[seq(1, length(hr_xval), by = 2)])
#for easting, we have to double the resolution to account for tiles#
hr_x <- ncdim_def('x', 'm', 
          vals = hr_xval,
          create_dimvar = T)
uv_geos_utm <- nc_open(filename = uv_path)

E_flux_var_t <- ncvar_def('E_flux', units = 'mW/ m-2',
                             dim = list(hr_x, uv_geos_utm$dim$y, date_dim),
                             missval = NA_real_, prec = 'float',
                          longname = "Daily Energy Flux (mW m-2) for 2022")

E_flux_var_mean <- ncvar_def('E_flux', units = 'mW/ m-2',
                        dim = list(hr_x, uv_geos_utm$dim$y),
                        missval = NA_real_, prec = 'float', longname = "Mean 2022 Energy Flux (mW m-2)")
crs_var <- uv_geos_utm$var$crs 
nc_close(uv_geos_utm)
E_flux_nc <- nc_create(filename = paste(r'(C:\Users\fbutschek\Desktop\hycom\2022\)', 
                                        'E_Flux_revisedv2.nc', sep = "/"), vars = list(E_flux_var_t,
                                                                             crs_var), 
                       force_v4 = TRUE)
E_flux_mean_nc <- nc_create(filename = paste(r'(C:\Users\fbutschek\Desktop\hycom\2022\)', 
                                        'mean_E_Flux_mean_revisedv2.nc', sep = "/"), 
                            vars = list(E_flux_var_mean, crs_var), 
                       force_v4 = TRUE)
E_flux_nc_blocked <- nc_create(filename = paste(r'(C:\Users\fbutschek\Desktop\hycom\2022\)', 
                                                'E_Flux_near_revised_blockedv2.nc', sep = "/"), 
                               vars = list(E_flux_var_t, crs_var), 
                               force_v4 = TRUE)
E_flux_mean_nc_blocked <- nc_create(filename = paste(r'(C:\Users\fbutschek\Desktop\hycom\2022\)', 
                                             'mean_E_Flux-blocked_mean_revisedv2.nc', sep = "/"), 
                            vars = list(E_flux_var_mean, crs_var), 
                            force_v4 = TRUE)
source('04_data_analysis/02_functions/nc_metadata_assignment_helpers.R')
nc_assign_geo(E_flux_mean_nc_blocked, "EPSG:32618")
nc_def_axis(E_flux_mean_nc_blocked, crs_type = 'projected')
nc_assign_geo(E_flux_nc_blocked, "EPSG:32618")
nc_def_axis(E_flux_nc_blocked, crs_type = 'projected')
nc_assign_geo(E_flux_mean_nc, "EPSG:32618")
nc_assign_geo(E_flux_nc, "EPSG:32618")
nc_def_axis(E_flux_mean_nc, crs_type = 'projected')
nc_def_axis(E_flux_nc, crs_type = 'projected')


on.exit(c(nc_close(E_flux_mean_nc),nc_close(E_flux_nc), 
          nc_close(E_flux_mean_nc_blocked), nc_close(E_flux_nc_blocked)), add = T)
n_nan <- 0
n_zero <- 0
nx_vec <- as.integer(stringr::str_split_i(names(out_ls3), pattern = '_', i = 2))
ny_vec <- as.integer(stringr::str_split_i(names(out_ls3), pattern = '_', i = 4))
n_tile <- as.integer(stringr::str_split_i(names(out_ls3), pattern = '_', i = 5))


# write values to nc_file -------------------------------------------------


for (i in 1:length(out_ls3)){
  metadat <- get(names(out_ls3[i]), fft_env)
  
  ny <- ny_vec[[i]]
  nx <- (nx_vec[[i]]-1)*2+n_tile[[i]] #increases res
  stopifnot("Tile indexing mismatch"=ny==metadat$ny & 
              nx_vec[[i]] ==metadat$nx & 
              n_tile[[i]]==metadat$tile_num)
    
  E_vec <-out_ls3[[i]]$E_flux
  mean_E <- mean(E_vec, na.rm = T)
  E_vec_blocked <-out_ls3[[i]]$E_flux_blocked
  mean_E_blocked <- mean(E_vec_blocked, na.rm = T)
  if(is.na(mean_E)) {print(names(out_ls3)[i])
    n_nan <- n_nan+1} # 148 nan are probably edge cases, TBC
  if(!is.na(mean_E) & dplyr::near(mean_E, 0)) n_zero <- n_zero+1
  # 113 zeros, TBC
  # e_mat[ny,nx] <- mean_E
  # e_array[ny,nx,]  <- E_vec
  ncvar_put(E_flux_nc, varid = "E_flux", vals = E_vec,
            start = c(nx,ny, 1), count = c(1,1,365))
  ncvar_put(E_flux_mean_nc,varid = "E_flux", vals = mean_E,
            start = c(nx,ny), count = c(1,1))
  
  ncvar_put(E_flux_nc_blocked, varid = "E_flux", vals = E_vec_blocked,
  start = c(nx,ny, 1), count = c(1,1,365))
  ncvar_put(E_flux_mean_nc_blocked,varid = "E_flux", vals = mean_E_blocked,
            start = c(nx,ny), count = c(1,1))
}



# add axis attributes -----------------------------------------------------


ncatt_put(E_flux_mean_nc, "E_flux", "grid_mapping", "crs")
ncatt_put(E_flux_mean_nc, "E_flux", "coordinates", "x y time")
ncatt_put(E_flux_nc, "E_flux", "grid_mapping", "crs")
ncatt_put(E_flux_nc, "E_flux", "coordinates", "x y")
nc_close(E_flux_mean_nc);nc_close(E_flux_nc)

ncatt_put(E_flux_mean_nc_blocked, "E_flux", "grid_mapping", "crs")
ncatt_put(E_flux_mean_nc_blocked, "E_flux", "coordinates", "x y time")
ncatt_put(E_flux_nc_blocked, "E_flux", "grid_mapping", "crs")
ncatt_put(E_flux_nc_blocked, "E_flux", "coordinates", "x y")
nc_close(E_flux_mean_nc_blocked);nc_close(E_flux_nc_blocked)


# convert to tif using {raster}--------------------------------------------
library(raster)
mean_flux_ncin <- nc_open(E_flux_mean_nc$filename)
mean_flux <- ncvar_get(mean_flux_ncin, varid = "E_flux")
# cell registration is bottom left for the nc files, which needs to be remapped
# extent is expanded to account fo centre cell registration in geotiff

x_min <- hr_xval[1]
x_max <- tail(hr_xval,1)+diff(tail(hr_xval,2))
y_vals <- E_flux_mean_nc$dim$y$vals
y_min <- min(y_vals) + diff(y_vals[1:2])
y_max <- max(y_vals) 

mean_flux_raster <- raster(t(mean_flux), xmn = x_min, xmx = x_max, 
                           ymn = y_min, ymx = y_max, 
                           crs = CRS("+init=EPSG:32618"))

res(mean_flux_raster) == c(unique(diff((mean_flux_ncin$dim$x$vals))),
                          -unique(diff((mean_flux_ncin$dim$y$vals))))
mean_flux_log10 <- log10(mean_flux_raster)
names(mean_flux_raster) <- "Mean.2022.Energy.Flux.mW.per.m2"
names(mean_flux_log10) <- paste('log10',names(mean_flux), sep = '.')
writeRaster(mean_flux_raster, 
            "C:\\Users\\fbutschek\\Desktop\\hycom\\2022\\mean_flux_revisedV2.tif",
            overwrite =T)
writeRaster(mean_flux_log10, 
            "C:\\Users\\fbutschek\\Desktop\\hycom\\2022\\log10mean_flux_revisedV2.tif",
            overwrite =T)

nc_pth <- r'(C:\Users\fbutschek\Desktop\hycom\2022\\mean_E_Flux-blocked_mean_revisedv2.nc)'
blocked_mean_flux_nc_in <- nc_open(nc_pth)
E_flux_mean_nc_blocked <- ncvar_get(blocked_mean_flux_nc_in, varid = "E_flux")
blocked_mean_flux_raster <- raster(t(E_flux_mean_nc_blocked), 
                                   xmn = x_min, xmx = x_max, 
                                   ymn = y_min, ymx = y_max, 
                                   crs = CRS("+init=EPSG:32618"))  
names(blocked_mean_flux_raster) <- "Mean.2022.Energy.Flux.mW.per.m2"
writeRaster(blocked_mean_flux_raster, 
            "C:\\Users\\fbutschek\\Desktop\\hycom\\2022\\blocked_mean_flux_revisedV2.tif",
            overwrite =T)
blocked_log10mean_flux_raster<- log10(blocked_mean_flux_raster)
names(blocked_log10mean_flux_raster) <- paste('log10',names(mean_flux), sep = '.')

writeRaster(blocked_log10mean_flux_raster, 
            "C:\\Users\\fbutschek\\Desktop\\hycom\\2022\\blocked_log10mean_flux_revisedV2.tif",
            overwrite =T)


fraction_blocked <- (mean_flux_raster-blocked_mean_flux_raster)/mean_flux_raster
names(fraction_blocked) <- "Mean.Blocked.Energy.Flux.Fraction.2022"
writeRaster(fraction_blocked, 
            "C:\\Users\\fbutschek\\Desktop\\hycom\\2022\\fraction_blocked_north_revised.tif",
            overwrite =T)


