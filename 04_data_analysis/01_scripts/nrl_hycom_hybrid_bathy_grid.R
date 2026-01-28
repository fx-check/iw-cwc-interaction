
# Create spatial mask from NRL bathymetry to match HYCOM dowlnoad#

# Preliminaries -----------------------------------------------------------
# 1) Import data and get names 
library(ncdf4)
# Having tested, the gebco grid is not a great match for HYCOM data, let's use the  
# Naval Research Laboratory (NRL) Digital Bathymetry Data Base (DBDB2) instead:
# download under 'topography' from: https://www.hycom.org/data/glba0pt08/expt-90pt8
# import NRL GLB bathymetry that the HYCOM model is based on
fn <- "02_secondary_data/depth_GLBa0.08_09.nc"
nc <- nc_open(fn)

# 1) Discover names for lon/lat in your file
vars <- names(nc$var)
vars

# Find lon/lat variable names
lon_name <- grepv('long', vars, ignore.case = T)
lat_name <- grepv('lat', vars, ignore.case = T)
depth_name <- vars[vars %in% c("depth","Depth","bathy","bathymetry")][1]
stopifnot(!is.na(lon_name), !is.na(lat_name), !is.na(depth_name))

# 2) Read 2-D lon/lat and depth matrices (dimensions should match)
plon <- ncvar_get(nc, lon_name)  # matrix [nj, ni] or [ni, nj] depending on file
plat <- ncvar_get(nc, lat_name)
depth <- ncvar_get(nc, depth_name)

nc_close(nc)

# Reduce global bathymetry to case study bounding box ---------------------
# If longitudes are 0..360, convert to -180..180
if (max(plon, na.rm = TRUE) > 180) plon <- ifelse(plon > 180, plon - 360, plon)

# 3) Build the per-cell mask from the 2-D matrices, here we use our study area
lon_min <- -80; lon_max <- -76 # modify this if needed
lat_min <-  30; lat_max <-  33 # modify this to match case study BBox

mask <- (plon >= lon_min & plon <= lon_max) &
  (plat >= lat_min & plat <= lat_max)

any(mask)       # should be TRUE
which(mask, arr.ind = TRUE)[1:5, , drop = FALSE]  # quick peek

# 4) Get the compact mask==TRUE region for further work in our study area
jj <- which(rowSums(mask, na.rm = TRUE) > 0)
ii <- which(colSums(mask, na.rm = TRUE) > 0)
j0 <- min(jj); j1 <- max(jj)
i0 <- min(ii); i1 <- max(ii)

plon_sub  <- t(plon [j0:j1, i1:i0]) # longitude increasing to the right
plat_sub  <- t(plat [j0:j1, i1:i0]) # latitude decreasing down 
depth_sub <- t(depth[j0:j1, i1:i0])
mask_sub  <- mask [j0:j1, i0:i1]
# Quick check
plon_sub[1:6,1:6] # looks correct
plat_sub[1:6,1:6] # needs to be decreasing to match raster expectation
depth_sub[1:6,1:6] # correct placement of na values

# Build a raster object from data subset ----------------------------------
library(raster)
fn2 <- choose.files() # choose a HYCOM .nc file that's been extracted with the same
#BBox coordinates, the following operation just imports the first var and z-dim.
hycom_grid <- raster(fn2) # which is fine for getting the xy dims:
e1 <- extent(hycom_grid)

dim(depth_sub) <- c(x = 44, y = 51) #add unequivocal dimensions along x=lon, y=lat
nrl_grid <- raster(depth_sub, 
       xmn = e1[1], xmx=e1[2], ymn = e1[3], ymx=e1[4], #extent
       crs = crs(hycom_grid)) #crs
image(nrl_grid) #checks out well
#additional sanity checks on the orientation of the NRL data:
# image(raster(plat_sub, 
#              xmn = e1[1], xmx=e1[2], ymn = e1[3], ymx=e1[4], #extent
#              crs = crs(hycom_grid))
# )
# image(raster(plon_sub, 
#              xmn = e1[1], xmx=e1[2], ymn = e1[3], ymx=e1[4], #extent
#              crs = crs(hycom_grid))
# )


# Resample to HYCOM resolution --------------------------------------------
# possibly due to the use of regional bathys in some areas of the HYCOM domain
# the resolution of the NRL grid resolution is slightly lower than the HYCOM data
# Let's resample to match the two
h_nrl_res <- resample(nrl_grid, hycom_grid, method = "ngb")

# Get max hycom depth level using stars -----------------------------------

library(stars)
hycom_array <- read_stars(fn2)
hycom_array <- units::drop_units(hycom_array)

#create depth dimension, :
hycom_array$depth <- st_dim_to_attr(hycom_array, which = which(names(st_dimensions(hycom_array))=='depth'))
# set depth values to NA when salinity is na:
hycom_array$depth[is.na(hycom_array$salinity)] <- NA
max_level <- apply(hycom_array$depth[,,,1],c(1,2),max, na.rm = T)
max_level[is.infinite(max_level)] <- NA #deal with infinite warning
max_level_raster <- raster(max_level,
                           xmn = e1[1], xmx=e1[2], ymn = e1[3], ymx=e1[4], #extent
                           crs = crs(hycom_grid)) #crs
depth_levels <- st_get_dimension_values(hycom_array, 
                                        which = which(names(st_dimensions(hycom_array))=='depth'))
depth_levels <- units::drop_units(depth_levels)
# Create hybrid raster form NRL and HYCOM max depth -----------------------

quantile(as.matrix(h_nrl_res)- t(max_level), p = seq(0,1,.1), na.rm = T)
# significant outliers on the extreme ends

# generate a H_true_hybrid array that is H_true or max_level + bin/2 when 
# max_level more than 1 depth bin shallower or deeper than H_true
H_hybrid <- h_nrl_res #initialize
# identify which bin is the deepest, by summing non-NA values
max_bin <- apply(hycom_array$depth,c(1,2), function(x){sum(!is.na(x))})
# if land, all NA --> will be zero

max_bin_grid <- raster(t(max_bin), xmn = e1[1], xmx=e1[2], ymn = e1[3], ymx=e1[4],
                  crs = crs(h_nrl_res))
image(max_bin_grid)
dz <- c(NA, diff(depth_levels)) #NA is for land here (0th bin)
rcl_mat <- data.frame('from'=cellStats(max_bin_grid, min, na.rm=T):
                        cellStats(max_bin_grid, max, na.rm=T),'to'=dz)
layer_width_below <- reclassify(max_bin_grid, rcl = rcl_mat) #layer width deeper than max cell
layer_width_above <- reclassify(max_bin_grid-1, rcl = rcl_mat) # layer width above max cell
#this works here because we don't have any cells in the 40th (5km) layer, otherwise need to add one
# entry to dz and rcl_mat

# Let's run some tests on differences:
cellStats(layer_width_below!=layer_width_above, stat = sum) # n of cells on the transition between layer widths

max_hycom_depth <- raster(t(max_level), xmn = e1[1], xmx=e1[2], ymn = e1[3], ymx=e1[4],crs = crs(H_true))
cellStats(h_nrl_res > max_hycom_depth + layer_width_below, sum)
cellStats(h_nrl_res < max_hycom_depth, sum) # quite a few cells shallower than HYCOM max depth
# add half the layer_width below to any NRL cells deepr than HYCOm max depth + full layer width
H_hybrid[h_nrl_res > max_hycom_depth + layer_width_below] <- 
  max_hycom_depth[h_nrl_res > max_hycom_depth + layer_width_below] + 
  layer_width_below[h_nrl_res > max_hycom_depth + layer_width_below] / 2
# for cells where gebco is too shallow, we'll leave the bathymetry as is and clip
# HYCOM data to true depth plus a 100m buffer (at maximum) to account for NRL
# uncertanties, as it is evident that the difference is < 50m in the vast majority:
hist(as.vector(H_hybrid[H_hybrid< max_hycom_depth])-
       as.vector(max_hycom_depth[H_hybrid< max_hycom_depth]))
# add 100m to max_hycom_depth to account for uncertainty in HYCOM and NRL data
# when NRL is too shallow
large_outliers <- (H_hybrid<= max_hycom_depth-100)
# cellStats(large_outliers, sum) #only 6 large outliers
small_outliers <- (H_hybrid> max_hycom_depth-100 & H_hybrid< max_hycom_depth)
H_hybrid[large_outliers] <- #the selection takes all the large outliers
  h_nrl_res[large_outliers]+100 #and accounts for NRL error
# The remaining data is much tighter around the zero difference:
# Here, just using the HYCOM max depth is appropriate
hist(as.vector(H_hybrid[small_outliers])-
       as.vector(max_hycom_depth[small_outliers]))
H_hybrid[small_outliers] <- max_hycom_depth[small_outliers]
# a small number of outliers will have shifted from >100 to <100 now, which is fine
cellStats(H_hybrid < max_hycom_depth,sum)
# The number of bathy data too shallow outlier cells has reduced by 700 to 6, 

# Next up, creating the mask using a depth weighted mean on a window function 






