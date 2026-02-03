# open first net_cdf file as read and write in the designated folder, add all other net cdf files as additional time steps
#
nc_dir <- r'(C:\Users\fbutschek\Desktop\hycom\2022\u_v_bottom)'
nc_files <- list.files(nc_dir, pattern = '.nc$', full.names = T)
library(ncdf4)
# create master file with structure 'like' nc_files[1]:
master_file <-nc_open(nc_files[1], write = F)
# redefine time dimension to accomodate for expaning this without limit:
time_dim <- ncdim_def(name = 'time', units = 'hours since 2000-01-01 00:00:00', vals = 1:length(nc_files), unlim = F)
# redefine variables to accomodate for expaning time dimension:
water_u_bottom_var <- ncvar_def(name = 'water_u_bottom', units = 'm/s', 
                                dim = list(ncdim_def('lon', 'degrees_east', vals = ncvar_get(master_file, 'lon')),
                                           ncdim_def('lat', 'degrees_north', vals = ncvar_get(master_file, 'lat')),
                                           time_dim), missval = -9999, prec = 'float')
water_v_bottom_var <- ncvar_def(name = 'water_v_bottom', units = 'm/s', 
                                dim = list(ncdim_def('lon', 'degrees_east', vals = ncvar_get(master_file, 'lon')),
                                           ncdim_def('lat', 'degrees_north', vals = ncvar_get(master_file, 'lat')),
                                           time_dim), missval = -9999, prec = 'float')


# create new master file with updated dimensions and variables:
master_file <- nc_create(filename = r'(C:\Users\fbutschek\Desktop\hycom\2022\u_v_bottom\master_file.nc)', 
                         vars = list(water_u_bottom_var, water_v_bottom_var)
                         )

                                
for (i in 1:length(nc_files)) {
  cat('Adding file', i, 'of', length(nc_files), '\n')
  current_file <- nc_open(nc_files[i])
  # get time variable to add as a dimension:
  time_var <- ncvar_get(current_file, 'time')
  # get two variables: 'water_u_bottom' and 'water_v_bottom':
  water_u_bottom <- ncvar_get(current_file, 'water_u_bottom')
  water_v_bottom <- ncvar_get(current_file, 'water_v_bottom')
  # append time var as dimension:
  ncvar_put(master_file, 'time', time_var, start = i, count = c(length(time_var)))
  # append  'water_u_bottom' and 'water_v_bottom' variables:
  ncvar_put(master_file, varid = 'water_u_bottom',
            start = c(1, 1, i), 
            vals = water_u_bottom,
            count = c(dim(water_u_bottom)[1], dim(water_u_bottom)[2], length(time_var))
            )
  ncvar_put(master_file, varid = 'water_v_bottom',
            start = c(1, 1, i), 
            vals = water_v_bottom,
            count = c(dim(water_u_bottom)[1], dim(water_u_bottom)[2], length(time_var))
  )
  nc_close(current_file)
}

nc_close(master_file)
uv_2022 <- read_stars(r'(C:\Users\fbutschek\Desktop\hycom\2022\u_v_bottom\master_file.nc)')
uv_2022
(daily_uv_2022 <- aggregate(uv_2022, '1 day', mean, na.rm = T))

(daily_uv_geos_2022 <-st_apply(daily_uv_2022,'time',
                              function(x) {
                                raster::movingFun(
                                  x,
                                  n = 5, #5 day interval at 24h timesteps 
                                  fun = mean,
                                  na.rm = TRUE
                                )
                              }) )
write_stars(daily_uv_geos_2022, dsn = paste(r'(C:\Users\fbutschek\Desktop\hycom\)',
                                            'daily_u_geos_2022.nc',sep = '/'), layer = 'water_u_bottom')

write_stars(daily_uv_geos_2022, dsn = paste(r'(C:\Users\fbutschek\Desktop\hycom\)',
                                            'daily_v_geos_2022.nc',sep = '/'), layer = 'water_v_bottom')
getwd()
