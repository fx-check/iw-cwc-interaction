
# Script to compile annual time series of hycom bottom u & v from 3-hourly .nc files

library(ncdf4)
nc_dir <- r'(C:\Users\fbutschek\Desktop\hycom\2022\u_v_bottom)'
nc_files <- list.files(nc_dir, pattern = '[0-9].nc$', full.names = T)

# create master file with structure 'like' nc_files[1]:
template_file <-nc_open(nc_files[5], write = F)
# redefine time dimension to accomodate for expanding this one-year limit:
time_dim <- ncdim_def(name = 'time', units = 'hours since 2000-01-01 00:00:00', 
                      vals = 1:length(nc_files), unlim = F)
#the actual time values are grabbed from the individual .nc files, this is just
# a place holder

# redefine variables to accomodate for expaning time dimension:
water_u_bottom_var <- ncvar_def(name = 'water_u_bottom',units = 'm/s', 
                                dim = list(ncdim_def('lon', 'degrees_east', 
                                                     vals = ncvar_get(template_file, 'lon')),
                                           ncdim_def('lat', 'degrees_north', 
                                                     vals = ncvar_get(template_file, 'lat')),
                                                     time_dim), 
                                missval = -9999, 
                                prec = 'float')
water_v_bottom_var <- ncvar_def(name = 'water_v_bottom', units = 'm/s', 
                                dim = list(ncdim_def('lon', 'degrees_east', 
                                                     vals = ncvar_get(template_file, 'lon')),
                                           ncdim_def('lat', 'degrees_north', 
                                                     vals = ncvar_get(template_file, 'lat')),
                                           time_dim), 
                                missval = -9999, 
                                prec = 'float')

# Define a'crs' variable:
crs_var <- ncvar_def("crs", "", list(), prec = "integer")


uv_geos <- nc_create(filename = r'(C:\Users\fbutschek\Desktop\hycom\2022\u_v_bottom\master_file.nc)', 
                         vars = list(water_u_bottom_var, water_v_bottom_var,
                                     crs_var),
                         force_v4 = T)
on.exit(nc_close(uv_geos), add = T)
glob_atts <- tail(ncatt_get(template_file, 0),6)
for(i in 1:length(glob_atts)){
  ncatt_put(uv_geos, 0, 
            attname = names(glob_atts[i]), 
            attval = glob_atts[[i]])
}

# Attach CF projection info (latitude-longitude on WGS84)
ncatt_put(uv_geos, "crs", "grid_mapping_name", "latitude_longitude")
ncatt_put(uv_geos, "crs", "longitude_of_prime_meridian", 0.0)
ncatt_put(uv_geos, "crs", "semi_major_axis", 6378137.0)
ncatt_put(uv_geos, "crs", "inverse_flattening", 298.257223563)


# WKT1 CRS definition:
wkt1 <- paste0(
  'GEOGCS["WGS 84",',
  'DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],',
  'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
)
ncatt_put(uv_geos, "crs", "spatial_ref", wkt1)
ncatt_put(uv_geos, 0, "Conventions", "CF-1.7")

# Link data variables to the CRS
ncatt_put(uv_geos, "water_u_bottom", "grid_mapping", "crs")
ncatt_put(uv_geos, "water_v_bottom", "grid_mapping", "crs")

for (i in 1:length(nc_files) ) {
  cat('Adding file', i, 'of', length(nc_files), '\n')
  current_file <- nc_open(nc_files[i])
  # get time variable to add as a dimension:
  time_var <- ncvar_get(current_file, 'time')
  # get two variables: 'water_u_bottom' and 'water_v_bottom':
  water_u_bottom <- ncvar_get(current_file, 'water_u_bottom')
  water_v_bottom <- ncvar_get(current_file, 'water_v_bottom')
  # append time var as dimension:
  ncvar_put(uv_geos, 'time', time_var, start = i, count = c(length(time_var)))
  # append  'water_u_bottom' and 'water_v_bottom' variables:
  ncvar_put(uv_geos, varid = 'water_u_bottom',
            start = c(1, 1, i), 
            vals = water_u_bottom,
            count = c(dim(water_u_bottom)[1], dim(water_u_bottom)[2], length(time_var))
            )
  ncvar_put(uv_geos, varid = 'water_v_bottom',
            start = c(1, 1, i), 
            vals = water_v_bottom,
            count = c(dim(water_u_bottom)[1], dim(water_u_bottom)[2], length(time_var))
  )
  nc_close(current_file)
}
# If you get a footprint that is 'apparently' enlarged by a half-cell in a GIS, 
# ensure to load the entire data set instead of just one variable, which may load
# incorrectly and cache stats in an .aux.xml., the defined crs should also help

# Aggregate to daily files and remove ageostrophic signals ----------------

# to take daily means and smooth out ageostrophic flows by a moving average,
# the stars library is most efficient. Caution is required between array/matrix 
# and geospatial conventions, with rows written ascending in the former and
# descending in the latter
library(stars)
nc_close(uv_geos) #required before loading with stars
uv_2022 <- read_stars(r'(C:\Users\fbutschek\Desktop\hycom\2022\u_v_bottom\master_file.nc)')
uv_2022
# downsample temporal resolution to daily using mean aggregation:
(daily_uv_2022 <- aggregate(uv_2022, '1 day', mean, na.rm = T))

# remove ageostrophic flows using 5-day moving average filter:
uv_geos_2022 <- st_apply(daily_uv_2022, MARGIN = c('x','y'),
                         FUN = function(x) {zoo::rollapply(x,
                                        width = 5, #5 days running mean
                                        FUN = mean, na.rm = T,
                                        partial = 3, # min 3-day mean for edges
                                        fill = NA
                         )}
)


# Rectify dimensions and write to nc file ---------------------------------
# here, we first re-assign the time dim, with dates lost through roll-apply
# also, we need to  permute and transpose the dimensions to account for the 
# aforementioned mismatch between geospatial and array indexing

st_dimensions(uv_geos_2022) <- st_dimensions(daily_uv_2022)
daily_uv_geos_2022 <- aperm(uv_geos_2022,c(2,3,1))
dim(daily_uv_geos_2022)#for some reason, this doesn't translate into the array
nr <- dim(daily_uv_geos_2022$water_u_bottom)['y']
geos_u_array <- daily_uv_2022$water_u_bottom[,,nr:1]
geos_u_array <- aperm(geos_u_array, c(2,3,1))
image(geos_u_array[,,123]) # looks good now
geos_v_array <- daily_uv_2022$water_v_bottom[,,nr:1]
geos_v_array <- aperm(geos_v_array, c(2,3,1))
image(geos_v_array[,,321]) # looks good now
#the following should be fairly similar:
test_u <- ncvar_get(template_file, 'water_u_bottom')
image(test_u) #3am data from second day
image(geos_u_array[,,2]) # second daily file, partly made up by test_u

#instead of writing these files with stars, which exports each timestep as a separate variable band,
# let's create an nc_file, the only vabeat is that we have to flip the nrow dimension to write correctly
u_geos <- ncvar_def('u_geostrophic', units = 'ms-1', missval = NA_real_, prec = 'float',
                    longname = "eastward geostrophic bottom current velocity",
                    dim = list(ncdim_def('lon', 'degrees_east', vals = ncvar_get(template_file, 'lon')),
                               ncdim_def('lat', 'degrees_north', vals = ncvar_get(template_file, 'lat')),
                               ncdim_def('time', 'days', vals = 0:364, 
                                         longname = 'days from 01-01-2022')
                               )
                    )
v_geos <- ncvar_def('v_geostrophic', units = 'ms-1', missval = NA_real_, prec = 'float',
                    longname = "northward geostrophic bottom current velocity",
                    dim = list(ncdim_def('lon', 'degrees_east', vals = ncvar_get(template_file, 'lon')),
                               ncdim_def('lat', 'degrees_north', vals = ncvar_get(template_file, 'lat')),
                               ncdim_def('time', 'days', vals = 0:364, 
                                         longname = 'days from 01-01-2022'))
                    )
uv_geos <- nc_create(paste(r'(C:\Users\fbutschek\Desktop\hycom\)','daily_uv_geos_2022.nc',
                sep = '/'), vars = list(u_geos, v_geos, crs_var), force_v4 = T)
on.exit(nc_close(uv_geos), add = T)

ncvar_put(uv_geos, varid = 'u_geostrophic', vals = geos_u_array,
          start = c(1,1,1), count = unname(dim(geos_u_array)))
ncvar_put(uv_geos, varid = 'v_geostrophic', vals = geos_v_array,
          start = c(1,1,1), count = dim(geos_u_array))

# Attach CF projection info (latitude-longitude on WGS84)
ncatt_put(uv_geos, "crs", "grid_mapping_name", "latitude_longitude")
ncatt_put(uv_geos, "crs", "longitude_of_prime_meridian", 0.0)
ncatt_put(uv_geos, "crs", "semi_major_axis", 6378137.0)
ncatt_put(uv_geos, "crs", "inverse_flattening", 298.257223563)


ncatt_put(uv_geos, "crs", "spatial_ref", wkt1)
ncatt_put(uv_geos, 0, "Conventions", "CF-1.7")

# Link data variables to the CRS
ncatt_put(uv_geos, "u_geostrophic", "grid_mapping", "crs")
ncatt_put(uv_geos, "v_geostrophic", "grid_mapping", "crs")
nc_close(master_file); nc_close(template_file); nc_close(uv_geos)

# The below method also writes to file, but each day is an own layer with Band1:Band365
# write_stars(daily_uv_geos_2022, dsn = paste(r'(C:\Users\fbutschek\Desktop\hycom\)',
#                                             'daily_u_geos_2022.nc',sep = '/'), layer = 'water_u_bottom')
# 
# write_stars(daily_uv_geos_2022, dsn = paste(r'(C:\Users\fbutschek\Desktop\hycom\)',
#                                             'daily_v_geos_2022.nc',sep = '/'), layer = 'water_v_bottom')
