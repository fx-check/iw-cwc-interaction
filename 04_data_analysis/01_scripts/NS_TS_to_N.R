
# script to compute N2 from WOA S-T data 
# library(stars)
library(ncdf4)
library(oce)
library(stars)

# Load data and get metadata ----------------------------------------------

in_dir <- r'(03_secondary_data\NovaScotia\intermediary_files)'
nc_Annual_t <- nc_open(paste(in_dir,'NS_T_Annual_temp.nc', sep = '/'))
Annual_t <- ncvar_get(nc_Annual_t, 't_an')
nc_Annual_sal <- nc_open(paste(in_dir,'NS_s_Annual_sal.nc', sep = '/'))
Annual_sal <- ncvar_get(nc_Annual_sal, 's_an')

nx <- nc_Annual_t$dim[["lon"]]$len
ny <- nc_Annual_t$dim[["lat"]]$len
K  <- nc_Annual_t$dim[["depth"]]$len

longs <- ncvar_get(nc_Annual_t, "lon")   
lats  <- ncvar_get(nc_Annual_t, "lat")    

## Initialise and calculate absolute salinity  ---------------------------


# generate matrices of long and lat, in [nx, ny] format native to nc & stars
# this is transposed from the visually more intuitive [ny, nx] format, sensu
# “rows are latitudes (ny), columns are longitudes (nx)”
lon_mat <- matrix(rep(longs,  times = ny), 
                  nrow = nx, 
                  ncol = ny, 
                  byrow = FALSE)
range(lats)
mid_lat <- mean(lats) #given 4.5 degree range of latitude, mid value is fine for p:
depth_levels <- nc_Annual_t$dim[["depth"]]$vals
p_levels <- gsw_p_from_z(z = -depth_levels, latitude = mid_lat)
# calculate absolute salinity from practical salinity, 
# Build lon/lat matrices matching x × y
# longs: length nx (x dimension), lats: length ny (y dimension)
lon_mat <- matrix(longs, nrow = ny, ncol = nx, byrow = TRUE)   # each column is 'longs'
lat_mat <- matrix(lats,  nrow = ny, ncol = nx, byrow = FALSE)  # each row has its own lat
# the above matrix computation is not strictly needed because gsw does it on the fly,
# but it will save computation time, also convert to positive [0,360]
lon_mat <- ifelse(lon_mat < 0, lon_mat + 360, lon_mat)
#initialize/pre-allocate array:
SA_arr <- array(NA_real_, dim = dim(Annual_sal))  # x × y × pressure
# create 3d arrays of long, lat and pressure
p3   <- array(rep(p_levels, each = nx*ny), dim = dim(Annual_sal))
lon3 <- array(lon_mat, dim = dim(Annual_sal))
lat3 <- array(lat_mat, dim = dim(Annual_sal))
SA_arr <- gsw::gsw_SA_from_SP(Annual_sal, p3, lon3, lat3)
dim(SA_arr) <- dim(Annual_sal)
chunks <- c(min(nx, 128), 
            min(ny, 128), 
            min(K, 5))
abs_salinity <- ncvar_def(name = 'absolute_salinity', units = 'g/kg',
                          dim = list(
                            nc_Annual_sal$dim$lon,
                            nc_Annual_sal$dim$lat,
                            ncdim_def('pressure', 'dbar', vals = p_levels)),
                          missval = NA, prec = 'float',
                          compression = 5L,
                          shuffle = F,
                          chunksizes = chunks,
                          longname    = "Sea water Absolute Salinity")


#store as a nc file:
abs_sa_nc <- nc_create(paste(in_dir,"Annual_abs_salinity.nc", sep ='/'), 
                       vars = abs_salinity, force_v4 = TRUE ) 
on.exit(nc_close(abs_sa_nc), add=TRUE)
ncvar_put(abs_sa_nc, "absolute_salinity",
          vals  = SA_arr,
          start = c(1, 1, 1),
          count = c(nx, ny, K))
nc_close(abs_sa_nc)

# repeat this framework for Conservative temperature:
# Conservative temperature from pressure and absolute salinity ------------
# Pre-allocate nc file with appropriate variable definition:
# create an nc file for Conservative temperature, using nc_Annual_t:
cons_temp <- ncvar_def(name = 'conservative_temperature', units = 'degC',
                       dim = list(
                         nc_Annual_t$dim$lon,
                         nc_Annual_t$dim$lat,
                         ncdim_def('pressure', 'dbar', vals = p_levels)),
                       missval = NA, prec = 'float',
                       compression = 5L,
                       shuffle = F,
                       chunksizes = chunks,
                       longname    = "Sea water Conservative Temperature")
ct_nc <- nc_create(paste(in_dir,"Annual_conservative_temperature.nc", sep ='/'), 
                   vars = cons_temp, force_v4 = TRUE) 
on.exit(nc_close(ct_nc), add=TRUE)
CT <- gsw::gsw_CT_from_t(SA_arr, Annual_t, p3)
dim(CT) <- dim(Annual_t)
ncvar_put(ct_nc, "conservative_temperature",
          vals  = CT,
          start = c(1, 1, 1),
          count = c(nx, ny, K)
          )

nc_close(ct_nc)

#rinse repeat the same for winter and summer (if desired), just swap file names


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

# sa_nc <- nc_open("02_secondary_data/NovaScotia/annual_abs_salinity.nc")
# ct_nc <- nc_open("02_secondary_data/NovaScotia/annual_conservative_temperature.nc")
# SA_arr <- ncvar_get(sa_nc, 'absolute_salinity')
# CT_arr <- ncvar_get(ct_nc, 'conservative_temperature')
# nc_close(ct_nc);nc_close(sa_nc)
mid_p <- p_levels[1:(length(p_levels)-1)] + diff(p_levels)/2
# initialise N2 nc_file:

n2 <- ncvar_def(name = 'N2', units = 's-2',
                       dim = list(
                         ncdim_def("lon", "degrees_east",  vals = longs, 
                                   create_dimvar = TRUE),
                         ncdim_def("lat", "degrees_north", vals = lats, 
                                   create_dimvar = TRUE),
                         ncdim_def('pressure', 'dbar', vals = mid_p,
                                   create_dimvar = T)),
                missval = NA, prec = 'float',
                compression = 5L,
                shuffle = F,
                chunksizes = chunks,
                longname    = "Squared buoyancy frequency N2")

# create N2:
N2_array <- array(NA_real_, dim = c(nx,ny, K-1))  # x × y × pressure 
n2_nc <- nc_create(paste(in_dir,"N2.nc", sep ='/'), 
                  vars = n2, force_v4 = TRUE)
on.exit(nc_close(n2_nc), add=TRUE)
# calculate N2 for each cell across the pressure levels:
for (n in seq_len(nx)) {
  for (m in seq_len(ny)) {
    lat <- lats[m]
    SA_col <- SA_arr[n, m, ]
    CT_col <- CT_arr[n, m, ]
    N2_array[n, m, ] <- compute_N2_profile(SA_col, CT_col, p_levels, lat)
  }
}
ncvar_put(n2_nc, "N2",
          vals  = N2_array,
          start = c(1, 1, 1),
          count = c(nx, ny, K-1))
#clean up:
nc_close(n2_nc)
rm(SA_arr)
rm(CT)

## Repeat for summer N2
#'summer' vars have no prefix here
sa_nc <- nc_open(paste(in_dir,"abs_salinity.nc", sep = '/'))
ct_nc <- nc_open(paste(in_dir,"conservative_temperature.nc", sep = '/'))
SA_arr <- ncvar_get(sa_nc, 'absolute_salinity')
CT_arr <- ncvar_get(ct_nc, 'conservative_temperature')
nc_close(ct_nc);nc_close(sa_nc)

N2_array <- array(NA_real_, dim = c(nx,ny, K-1)) 
n2_nc <- nc_create(paste(in_dir,"summerN2.nc", sep ='/'), 
                   vars = n2, force_v4 = TRUE)
on.exit(nc_close(n2_nc), add=TRUE)
# calculate N2 for each cell across the pressure levels:
for (n in seq_len(nx)) {
  for (m in seq_len(ny)) {
    lat <- lats[m]
    SA_col <- SA_arr[n, m, ]
    CT_col <- CT_arr[n, m, ]
    N2_array[n, m, ] <- compute_N2_profile(SA_col, CT_col, p_levels, lat)
  }
}
ncvar_put(n2_nc, "N2",
          vals  = N2_array,
          start = c(1, 1, 1),
          count = c(nx, ny, K-1))
#clean up:
nc_close(n2_nc)
rm(SA_arr)
rm(CT_arr)
## Repeat for winter N2
sa_nc <- nc_open(paste(in_dir,"winter_abs_salinity.nc", sep = '/'))
ct_nc <- nc_open(paste(in_dir,"winter_conservative_temperature.nc", sep = '/'))
SA_arr <- ncvar_get(sa_nc, 'absolute_salinity')
CT_arr <- ncvar_get(ct_nc, 'conservative_temperature')
nc_close(ct_nc);nc_close(sa_nc)
N2_array <- array(NA_real_, dim = c(nx,ny, K-1)) 
n2_nc <- nc_create(paste(in_dir,"winterN2.nc", sep ='/'), 
                   vars = n2, force_v4 = TRUE)
on.exit(nc_close(n2_nc), add=TRUE)
# calculate N2 for each cell across the pressure levels:
for (n in seq_len(nx)) {
  for (m in seq_len(ny)) {
    lat <- lats[m]
    SA_col <- SA_arr[n, m, ]
    CT_col <- CT_arr[n, m, ]
    N2_array[n, m, ] <- compute_N2_profile(SA_col, CT_col, p_levels, lat)
  }
}
ncvar_put(n2_nc, "N2",
          vals  = N2_array,
          start = c(1, 1, 1),
          count = c(nx, ny, K-1))
#clean up:
nc_close(n2_nc)
rm(SA_arr)
rm(CT)

# Max depth raster using stars [optional] ---------------------------------

library(stars)
woa_array <- read_stars(n2_nc$filename)
woa_array <- units::drop_units(woa_array)
woa_array <- st_redimension(woa_array, new_dims = st_dimensions(woa_array)[1:3])
#create depth dimension, :
# st_dim
woa_array$pressure <- st_dim_to_attr(woa_array, which = which(names(st_dimensions(woa_array))=='pressure'))
# set pressure values to NA when salinity is na:
woa_array$pressure[is.na(woa_array$N2.nc)] <- NA
max_level <- apply(woa_array$pressure,c(1,2),max, na.rm = T)
max_level[!is.finite(max_level)] <- NA #deal with infinite warning
max_level <- max_level[,ncol(max_level):1] #flipping required for nc and matrix
N2_levels <- mid_p
mid_levels <- p_levels[-c(1,K)]
rm(woa_array) #stars array in memory slow down the system
# Build bottom weights and bottom N2 --------------------------------------
source("04_data_analysis/02_functions/bottom_window_weights_N2_no_grid.R")
bottom_weights <- build_n2_bottom_weights(max_n2_depth = max_level, depth_levels = mid_p)
bottom_n2 <- ncvar_def(name = 'Bottom N2', units = 's^-2',
                       dim = list(
                         nc_Annual_t$dim$lon,
                         nc_Annual_t$dim$lat),
                       missval = NA, prec = 'float',
                       compression = 5L,
                       shuffle = F,
                       longname    = "Squared Annual bottom buoyancy frequency N2")
bottom_n2_nc <- nc_create(paste(in_dir,"Annual bottom_N2.nc", sep ='/'), 
                          vars = bottom_n2, force_v4 = TRUE )
on.exit(nc_close(bottom_n2_nc), add=TRUE)

# sanity check:
apply(bottom_weights[,,],MARGIN = c(1,2), sum) #looks good
all.equal(dim(bottom_weights), dim(N2_array)) #sanity check
bottom_weights[which(is.na(N2_array[N2_array[,,1]]))]
weighted_N2_array <- N2_array* bottom_weights
bottom_n2_matrix <- apply(weighted_N2_array,MARGIN = c(1,2), sum, na.rm = T)
ncvar_put(bottom_n2_nc, "Bottom N2",
          vals  = bottom_n2_matrix,
          start = c(1, 1),
          count = c(nx, ny))
nc_close(bottom_n2_nc)
