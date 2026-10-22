
# script to compute N2 from high-dimensional modelled S-T data 
library(stars)
library(ncdf4)
st_2022 <- read_stars(r'(C:\Users\fbutschek\Desktop\hycom\2022_south\st_master_file.nc)',
                      proxy = T)
library(oce)

# Initialise and calculate absolute salinity ------------------------------
lats <- st_get_dimension_values(st_2022, 
                                which(names(st_dimensions(st_2022))=='y'))
mid_lat <- mean(lats) #given 3 degree range of latitude, mid value is fine
longs <- st_2022 |> # get and convert longs to + [0,360] E of Greenwich
  st_get_dimension_values(,which = which(names(dim(st_2022))=='x')
  )%% 360 
depth_levels <- st_get_dimension_values(st_2022, 
                                        which(names(st_dimensions(st_2022))=='depth'))
depth_levels <- units::drop_units(depth_levels)
p_levels <- gsw_p_from_z(z = -depth_levels, latitude = mid_lat)
st_2022 <- st_set_dimensions(st_2022, 'depth', 
                             values = p_levels) #change depth to pressure levels
names(st_dimensions(st_2022) #rename depth to pressure:
      )[which(names(st_dimensions(st_2022))=='depth')] <- 'pressure' 
dims_xykt <- dim(st_2022)
nx <- dims_xykt[1]
ny <- dims_xykt[2]
K  <- dims_xykt[3]
tm <- dims_xykt[4]
# calculate absolute salinity from practical salinity, 
# Build lon/lat matrices matching x × y
# longs: length nx (x dimension), lats: length ny (y dimension)
lon_mat <- matrix(longs, nrow = ny, ncol = nx, byrow = TRUE)   # each column is 'longs'
lat_mat <- matrix(lats,  nrow = ny, ncol = nx, byrow = FALSE)  # each row has its own lat
# the above matrix computation is not strictly needed because gsw does it on the fly,
# but it will save computation time

# This has already been done in the calculation of 'longs' above, but sanity check:
lon_mat <- ifelse(lon_mat < 0, lon_mat + 360, lon_mat)
# Pre-allocate
SA_arr <- array(NA_real_, dim = dim(st_2022))  # x × y × pressure × time

# compute time-step-wise absolute salinity --------------------------------
nc <- nc_open(r'(C:\Users\fbutschek\Desktop\hycom\2022_south\st_master_file.nc)')
on.exit(nc_close(nc), add=TRUE)

# sizes
nx <- nc$dim[["lon"]]$len
ny <- nc$dim[["lat"]]$len
K  <- nc$dim[["depth"]]$len
tm <- nc$dim[["time"]]$len

longs <- ncvar_get(nc, "lon")   
lats  <- ncvar_get(nc, "lat")    
# generate matrices of long and lat, in [nx, ny] format native to nc & stars
# this is transposed from the visually more intuitive [ny, nx] format, sensu
# “rows are latitudes (ny), columns are longitudes (nx)”
lon_mat <- matrix(rep(longs,  times = ny), 
                  nrow = nx, 
                  ncol = ny, 
                  byrow = FALSE)
lat_mat <- matrix(rep(lats,   each  = nx), nrow = nx, ncol = ny, byrow = FALSE)

sal_str <-  ncvar_get(nc, "salinity",
                      start = c(1, 1, 1, 1),       # x, y, pressure, time
                      count = c(nx, ny, K, 1)) #structure of time slice [nx, ny, K]
# create 3d arrays of long, lat and pressure
p3   <- array(rep(p_levels, each = nx*ny), dim = dim(sal_str))
lon3 <- array(lon_mat, dim = dim(sal_str))
lat3 <- array(lat_mat, dim = dim(sal_str))
#create absolute salinity variable definition:
#for comprssion:
chunks <- c(min(nx, 128), 
            min(ny, 128), 
            min(K, 5), 1) #chunksize is one time step

abs_salinity <- ncvar_def(name = 'absolute_salinity', units = 'g/kg',
                          dim = list(
                            nc$dim$lon,
                            nc$dim$lat,
                            ncdim_def('pressure', 'dbar', vals = p_levels),
                            nc$dim$time),
                          missval = NA, prec = 'float',
                          compression = 5L,
                          shuffle = F,
                          chunksizes = chunks,
                          longname    = "Sea water Absolute Salinity")

#store as a nc file:
abs_sa_nc <- nc_create(r'(C:\Users\fbutschek\Desktop\hycom\2022_south\abs_salinity.nc)', 
                       vars = abs_salinity, force_v4 = TRUE ) 
on.exit(nc_close(abs_sa_nc), add=TRUE)

for (i in seq_len(tm)) {
  SP_i <- ncvar_get(nc, "salinity",
                    start = c(1, 1, 1, i),       # x, y, pressure, time
                    count = c(nx, ny, K, 1))
  
  # Compute absolute salinity for this time slice
  SA_i <- gsw::gsw_SA_from_SP(SP_i, p3, lon3, lat3)
  # Write just this time-slice and move on
  ncvar_put(abs_sa_nc, "absolute_salinity",
            vals  = SA_i,
            start = c(1, 1, 1, i),
            count = c(nx, ny, K, 1))
}
#assign SA_arr to abs_sa_nc:
nc_close(abs_sa_nc)
rm(SA_i)

# repeat this framework for Conservative temperature:
# Conservative temperature from pressure and absolute salinity ------------

# Pre-allocate nc file with appropriate variable definition:
# create an nc file for Conservative temperature, using nc as blueprint:
cons_temp <- ncvar_def(name = 'conservative_temperature', units = 'degC',
                          dim = list(
                            nc$dim$lon,
                            nc$dim$lat,
                            ncdim_def('pressure', 'dbar', vals = p_levels),
                            nc$dim$time),
                          missval = NA, prec = 'float',
                          compression = 5L,
                          shuffle = F,
                          chunksizes = chunks,
                          longname    = "Sea water Conservative Temperature")
ct_nc <- nc_create(r'(C:\Users\fbutschek\Desktop\hycom\2022_south\ct.nc)', 
                       vars = cons_temp, force_v4 = TRUE ) 
on.exit(nc_close(ct_nc), add=TRUE)
# here, we may have to iterate through time steps and pressure levels
sa_nc <- nc_open(filename = r'(C:\Users\fbutschek\Desktop\hycom\2022_south\abs_salinity.nc)')
on.exit(nc_close(sa_nc), add=TRUE)

for (i in seq_len(tm)) {
  SA_i <- ncvar_get(sa_nc, "absolute_salinity",
                    start = c(1, 1, 1, i),       # x, y, pressure, time
                    count = c(nx, ny, K, 1))
  temp_i <- ncvar_get(nc, "water_temp",
                      start = c(1, 1, 1, i),       # x, y, pressure, time
                      count = c(nx, ny, K, 1))
  # p_levels are stored along the same dimensions in p3
  CT_i <- gsw::gsw_CT_from_t(SA_i, temp_i, p3)
  dim(CT_i) <- dim(temp_i)
  
  ncvar_put(ct_nc, "conservative_temperature",
            vals  = CT_i,
            start = c(1, 1, 1, i),
            count = c(nx, ny, K, 1))
}
nc_close(ct_nc)
nc_close(sa_nc)
nc_close(nc)

# Helper function to compute N2 profiles ----------------------------------

compute_N2_profile <- function(SA, CT, p, lat) {
  stopifnot(length(SA) == length(CT), length(CT) == length(p))
  K <- length(p)
  out <- rep(NA_real_, K - 1)  # result aligned to original mid-indices 1..K-1
  
  ok <- is.finite(SA) & is.finite(CT) & is.finite(p)
  if (!any(ok)) return(out) #abort if all NAs
  
  # Find contiguous TRUE runs in 'ok', here we ensure that adjacent levels are valid
  r <- rle(ok) #run length encoding to prevent out-of-grid mid-point around spurious NAs
  ends   <- cumsum(r$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  # Take valid runs with min. 2 levels to form 1+ pairs
  good_runs <- which(r$values & (r$lengths > 1L))
  for (j in good_runs) {
    a <- starts[j]
    b <- ends[j]
    N2 <- gsw::gsw_Nsquared(SA[a:b], CT[a:b], p[a:b], lat = lat)
    # Assign by original positions: pairs are a:(b-1)
    out[a:(b - 1L)] <- N2$N2
    # We deliberately don't use N2$p_mid for assignment to avoid float/name issues
  }
  out
}
# N2 calculation ----------------------------------------------------------

sa_nc <- nc_open("C:\\Users\\fbutschek\\Desktop\\hycom\\2022_south\\abs_salinity.nc")
on.exit(nc_close(sa_nc), add=TRUE)
ct_nc <- nc_open("C:\\Users\\fbutschek\\Desktop\\hycom\\2022_south\\ct.nc")
on.exit(nc_close(ct_nc), add=TRUE)
mid_p <- p_levels[1:(length(p_levels)-1)] + diff(p_levels)/2
# initialise N2 nc_file:
n2 <- ncvar_def(name = 'N2', units = 's^-2',
                       dim = list(
                         nc$dim$lon,
                         nc$dim$lat,
                         ncdim_def('pressure', 'dbar', vals = mid_p),
                         nc$dim$time),
                       missval = NA, prec = 'float',
                       compression = 5L,
                       shuffle = F,
                       chunksizes = chunks,
                       longname    = "Squared buoyancy frequency N2")

# create N2:
n2_nc <- nc_create(r'(C:\Users\fbutschek\Desktop\hycom\2022_south\N2.nc)', 
                  vars = n2, force_v4 = TRUE )
on.exit(nc_close(n2_nc), add=TRUE)
# for every timestep, load an array of sa and ct with all pressure levels
# calculate N2 for each cell across the pressure levels:
for (i in seq_len(tm)) {
  sa_xyk <- ncvar_get(sa_nc, "absolute_salinity",
                      start = c(1, 1, 1, i),       # x, y, pressure, time
                      count = c(nx, ny, K, 1))
  ct_xyk <- ncvar_get(ct_nc, "conservative_temperature",
                      start = c(1, 1, 1, i),       # x, y, pressure, time
                      count = c(nx, ny, K, 1))
  N2_array <- array(NA_real_, dim = c(nx, ny, K - 1L))
  
  for (n in seq_len(nx)) {
    for (m in seq_len(ny)) {
      lat <- lats[m]
      
      SA_col <- sa_xyk[n, m, ]
      CT_col <- ct_xyk[n, m, ]
      
      N2_array[n, m, ] <- compute_N2_profile(SA_col, CT_col, p_levels, lat)
    }
  }
  ncvar_put(n2_nc, "N2",
            vals  = N2_array,
            start = c(1, 1, 1, i),
            count = c(nx, ny, K-1, 1))
}


# Bottom N2 weights -------------------------------------------------------

# having computed the N2 raster for all timesteps and x,y cells, we next 
# calculate the bottom N2 weights based on the true depth raster and max depth levels
library(raster)
# source('04_data_analysis/02_functions/bottom_window_weights.R')

source('04_data_analysis/02_functions/bottom_window_weights_N2.R')

H_hybrid <- raster('02_secondary_data/NRL_HYCOM_hybrid_bathy_south.tif')
max_hycom_depth <- raster('03_derived_data/HYCOM_max_depth_south.tif')

dz1 <- diff(depth_levels) #computes inter depth level differences
# here we assume 0 top edge (surface) and bottom edge is last mid + half last diff:
depth_edges <- c(0, head(depth_levels, -1) + dz1/2, tail(depth_levels, 1) + tail(dz1, 1)/2)
dz2 <- diff(depth_edges)
# generate max N2 valid depth raster from nc file:
n2_ex <- brick(r'(C:\Users\fbutschek\Desktop\hycom\2022_south\N2.nc)', lvar = 4)
n2_max_valid <- calc(!is.na(n2_ex), fun = sum) #unambgiuous index of max N2 depth level 
n2_max_depth <- reclassify(n2_max_valid, 
                  rcl = cbind('is'=0:39, #we never get beyond 38, as 4000m the deepest cell in study area
                        'becomes'=  c(NA, # 0 becomes NA
                                      
                           depth_levels[1:39]+dz1[1:39]/2))) #deeper becomes mid value for N2
mid_levels <- depth_levels[1:39]+dz1[1:39]/2
# mid_levels <-  depth_levels[1:38]+diff(depth_levels[-length(depth_levels)])/2
hybrid_nc <- nc_open(filename = '02_secondary_data/NRL_HYCOM_hybrid_bathy_south.nc')
on.exit(nc_close(hybrid_nc), add = T)
names(hybrid_nc$var)
hybrid_depths <- ncvar_get(hybrid_nc, names(hybrid_nc$var)[2]
                           # 'NRL_HYCOM_hybrid_bathy'
                           )
n2_max_depth <- t(as.matrix(n2_max_depth))
all.equal(is.na(n2_max_depth), is.na(hybrid_depths))
# Any mismatches?
sum(is.na(t(as.matrix(n2_max_depth))))
#N:38, S: 0
sum(is.na(hybrid_depths))
#N: 48, S:0
n2_max_depth[is.na(hybrid_depths)] # an additional na 10 cells in very shallow water

bottom_weights <- build_n2_bottom_weights(hybrid_depths, n2_max_depth, depth_levels = mid_levels)
K <- length(mid_levels)
# bottom_weights <- aperm(bottom_weights,c(2,1,3))
# dim(bottom_weights)
# the generated bottom weights mask has zdim with 41 cells, i.e. two more than any N2 data
#
#initialize a nx × ny × K = 1  × time nc  raster:
bottom_n2 <- ncvar_def(name = 'Bottom N2', units = 's^-2',
                       dim = list(
                         nc$dim$lon,
                         nc$dim$lat,
                         nc$dim$time),
                       missval = NA, prec = 'float',
                       compression = 5L,
                       shuffle = F,
                       longname    = "Squared bottom buoyancy frequency N2")
bottom_n2_nc <- nc_create(r'(C:\Users\fbutschek\Desktop\hycom\2022_south\bottom_N2.nc)', 
                   vars = bottom_n2, force_v4 = TRUE )
on.exit(nc_close(bottom_n2_nc), add=TRUE)

# loop through time and calculate bottom N2 from n2_nc:
# sanity check:
apply(bottom_weights[,,],MARGIN = c(1,2), sum) #looks good

#the following loop really should be parallelised
for (i in seq_len(tm)) {
  n2_xyk <- ncvar_get(n2_nc, "N2",
                      start = c(1, 1, 1, i),       # x, y, pressure, time
                      count = c(nx, ny, K, 1))
  weighted_N2_array <- n2_xyk * bottom_weights[,ny:1,1:K] #reverse y-dim on purpose to align
  bottom_n2_matrix <- apply(weighted_N2_array,MARGIN = c(1,2), sum, na.rm = T)

  ncvar_put(bottom_n2_nc, "Bottom N2",
            vals  = bottom_n2_matrix,
            start = c(1, 1, i),
            count = c(nx, ny, 1))
}

bottom_n2_matrix <-ncvar_get(bottom_n2_nc, "Bottom N2",
                             start = c(1, 1, 50),
                             count = c(nx, ny,1))
nc_close(bottom_n2_nc)
image(raster(t(bottom_n2_matrix), 
             xmn = e1[1], xmx=e1[2], ymn = e1[3], ymx=e1[4], #extent
             crs = crs(hycom_grid)))


# # generate seasonal average using ncdf4 : --------------------------------
# the stars method didn't work
library(lubridate)
bottom_n2 <- nc_open(r'(C:\Users\fbutschek\Desktop\hycom\2022_south\bottom_N2.nc)')
time_seconds <- ncvar_get(bottom_n2, varid = 'time')*3600
time_origin <-  as_datetime("2000-01-01 00:00:00", tz= "UTC")
add_time <- as.numeric(24*3600*(time_origin-origin))
t_var <- lubridate::as_datetime(add_time + time_seconds)
(spring_start <- which(t_var == "2022-03-22 00:00:00 UTC"))
# last day of summer:
summer_end <- tail(which(as_date(t_var)=="2022-09-21"),1)
#filter time to before March 22 or after Sept 21:
late_winter_n2 <- ncvar_get(bottom_n2, 
                             start = c(1,1,1),
                             count = c(nx,ny,spring_start-1))
mean_lw_n2 <- apply(late_winter_n2, MARGIN = c(1,2), FUN = mean, na.rm = T)

summer_n2 <- ncvar_get(bottom_n2, 
                       start = c(1,1,spring_start),
                       count = c(nx,ny,summer_end-spring_start))
mean_summer_n2 <- apply(summer_n2, MARGIN = c(1,2), FUN = mean, na.rm = T)

early_winter_n2 <- ncvar_get(bottom_n2, 
                             start = c(1,1,summer_end+1),
                             count = c(nx,ny,length(t_var) - summer_end))
mean_ew_n2 <- apply(early_winter_n2, MARGIN = c(1,2), FUN = mean, na.rm = T)

mean_winter_n2 <- (mean_ew_n2+mean_lw_n2)/2


# Write both the mean_summer_n2 and mean_winter_n2 to a new .nc file as separate bands:
winter_n2nc <- ncvar_def(name = 'winter_N2', units = 's^-2',
                       dim = list(
                         nc$dim$lon,
                         nc$dim$lat),
                       missval = NA, prec = 'float',
                       compression = 5L,
                       shuffle = F#,
                       # longname    = "Winter (before 2022-03-22 and after 2022-09-21) average of squared bottom buoyancy frequency N2"
                       )
summer_n2nc <- ncvar_def(name = 'summer_N2', units = 's^-2',
                         dim = list(
                           nc$dim$lon,
                           nc$dim$lat),
                         missval = NA, prec = 'float',
                         compression = 5L,
                         shuffle = F#,
                         # longname    = "Summer (2022-03-22 to 2022-09-21) average of squared bottom buoyancy frequency N2"
                         )


seasonal_n2_nc <- nc_create(r'(C:\Users\fbutschek\Desktop\hycom\2022_south\seasonal_bottom_N2_data.nc)',
                            vars = list(winter_n2nc, summer_n2nc))
on.exit(nc_close(seasonal_n2_nc), add = T)
ncvar_put(seasonal_n2_nc, varid = 'winter_N2',
          vals =mean_winter_n2)
ncvar_put(seasonal_n2_nc, varid = 'summer_N2',
          vals = mean_summer_n2)
nc_close(seasonal_n2_nc)


# Reproject N2 bottom layer data to UTM18N --------------------------------
library(ncdf4); library(stars)

nc_grid <- rast(seasonal_n2$filename)
nc_grid_utm <- terra::project(nc_grid, "EPSG:32618", 
                              res = c(7450, 4380),
                              method = 'bilinear')
## Define target wrap dimension from nc template grid ----------------------
target <- st_as_stars(raster(nc_grid_utm))


winter_n2_stars <- read_stars(seasonal_n2_nc$filename, sub = 1)
st_crs(winter_n2_stars) <- "EPSG:4326"
winter_n2_utm <- st_warp(winter_n2_stars, target,
                           method = 'bilinear', use_gdal = T,
                           no_data_value = -9999)
names(winter_n2_utm) <- 'winter_n2'
summer_n2_stars <- read_stars(seasonal_n2_nc$filename, sub = 2)
st_crs(summer_n2_stars) <- "EPSG:4326"
summer_n2_utm <- st_warp(summer_n2_stars, target,
                         method = 'bilinear', use_gdal = T,
                         no_data_value = -9999)  
names(summer_n2_utm) <- 'summer_n2'

# write_stars(summer_n2_utm,r'(C:\Users\fbutschek\Desktop\hycom\2022_south\N2_UTMtest.tif)')
# Build NC files ----------------------------------------------------------

x_dim <- st_dimensions(target)$x
x_vals <- seq(from = x_dim$offset, by = x_dim$delta, length.out =x_dim$to)
# y_vals <- st_get_dimension_values(u_geos_2022_UTM18N, 'y')
y_dim <-st_dimensions(target)$y
y_vals <- seq(from = y_dim$offset, by = y_dim$delta, length.out =y_dim$to)

crs_var <- ncvar_def("crs", "", list(), prec = "integer")
source('04_data_analysis/02_functions/nc_metadata_assignment_helpers.R')

winter_n2_utm_nc <- ncvar_def(name = 'winter_N2', units = 's^-2', 
                           dim = list(ncdim_def('x', 'm', 
                                                vals = x_vals,
                                                create_dimvar = T),
                                      ncdim_def('y', 'm', 
                                                vals = y_vals,
                                                create_dimvar = T)),
                         missval = NA, prec = 'float',
                         compression = 5L,
                         shuffle = F,
                         longname    = "Winter (before 2022-03-22 and after 2022-09-21) average of squared bottom buoyancy frequency N2")
summer_n2_utm_nc <- ncvar_def(name = 'summer_N2', units = 's^-2', 
                         dim = list(ncdim_def('x', 'm', 
                                              vals = x_vals,
                                              create_dimvar = T),
                                    ncdim_def('y', 'm', 
                                              vals = y_vals,
                                              create_dimvar = T)), 
                         missval = NA, prec = 'float',
                         compression = 5L,
                         shuffle = F,
                         longname    = "Summer (2022-03-22 to 2022-09-21) average of squared bottom buoyancy frequency N2")

seasonal_n2_UTM <- nc_create(r'(C:\Users\fbutschek\Desktop\hycom\2022_south\seasonal_bottom_N2_data_UTM.nc)',
                            vars = list(winter_n2_utm_nc, summer_n2_utm_nc, crs_var))
nc_assign_geo(seasonal_n2_UTM, refsys = "EPSG:32618")
ncvar_put(seasonal_n2_UTM, 'winter_N2',
          vals = winter_n2_utm$winter_n2)
ncvar_put(seasonal_n2_UTM, 'summer_N2',
          vals = summer_n2_utm$summer_n2)

nc_def_axis(seasonal_n2_UTM, 'projected')

ncatt_put(seasonal_n2_UTM, "winter_N2", "grid_mapping", "crs")
ncatt_put(seasonal_n2_UTM, "summer_N2", "grid_mapping", "crs")

            
            
nc_close(seasonal_n2_UTM)
            