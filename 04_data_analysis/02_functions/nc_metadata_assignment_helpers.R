# Helper functions to reduce text for assigning geographic data to nc file
wkt_32618 <- paste0('PROJCS["WGS_1984_UTM_Zone_18N",',
               'GEOGCS["GCS_WGS_1984",',
                      'DATUM["D_WGS_1984",',
                          'SPHEROID["WGS_1984",6378137.0,298.257223563]],',
                      'PRIMEM["Greenwich",0.0],',
                      'UNIT["Degree",0.0174532925199433]],',#end of geo crs
               'PROJECTION["Transverse_Mercator"],',
               'PARAMETER["False_Easting",500000.0],',
               'PARAMETER["False_Northing",0.0],',
               'PARAMETER["Central_Meridian",-75.0],',
               'PARAMETER["Scale_Factor",0.9996],',
               'PARAMETER["Latitude_Of_Origin",0.0],',
               'UNIT["Meter",1.0]]')
  
# WKT1 CRS definition:
wkt_4326 <- paste0(
  'GEOGCS["WGS 84",',
  'DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],',
  'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]')
nc_assign_geo <- function(nc, refsys = 'EPSG:4326'){
  if(refsys == 'EPSG:32618'){
    ncatt_put(nc, "crs", "grid_mapping_name", "transverse_mercator")
    ncatt_put(nc, "crs", "false_easting", 500000.0)
    ncatt_put(nc, "crs", "false_northing", 0.0)
    ncatt_put(nc, "crs", "longitude_of_central_meridian", -75.0)
    ncatt_put(nc, "crs", "scale_factor_at_central_meridian", 0.9996)
    ncatt_put(nc, "crs", "latitude_of_projection_origin", 0.0)
    #for WGS84 geoid components
    ncatt_put(nc, "crs", "semi_major_axis", 6378137.0)
    ncatt_put(nc, "crs", "inverse_flattening", 298.257223563)
    #spatial ref:
    ncatt_put(nc, "crs", "spatial_ref", wkt_32618)
    ncatt_put(nc, "crs", "epsg_code", "EPSG:32618")
    
      } else if(refsys == 'EPSG:4326'){
    ncatt_put(nc, "crs", "grid_mapping_name", "latitude_longitude")
    ncatt_put(nc, "crs", "longitude_of_prime_meridian", 0.0)
    ncatt_put(nc, "crs", "semi_major_axis", 6378137.0)
    ncatt_put(nc, "crs", "inverse_flattening", 298.257223563)
    ncatt_put(nc, "crs", "spatial_ref", wkt_4326)
    ncatt_put(nc, "crs", "epsg_code", "EPSG:4326")
  }
  else {print('Selected refsys is not implemented, use "EPSG:4326" or "EPSG:32618" or recode function')}
}

nc_def_axis <- function(nc, crs_type = c('geographic', 'projected')){
  x_names <- c('x','X','lon', 'long', 'longitude','Lon', 'Long', 'Longitude')
  x_name <- x_names[x_names%in%names(nc$dim)]
  y_names <- c('y', 'Y','lat', 'latitude','Lat', 'Latitude')
  y_name <- y_names[y_names%in%names(nc$dim)]
  crs_type = crs_type[1]
  if(crs_type=='geographic'&& x_name %in% c('x','X') ||
     crs_type=='projected'&& x_name %in% 
     c('lon', 'long', 'longitude','Lon', 'Long', 'Longitude')){
    warning(paste0('Mismatch between nc dimnames and expected format, given selected crs type: ', 
                  crs_type, ', x_name:', x_name)
    )}
  if(crs_type=='geographic'&& y_name %in% c('y','Y') ||
     crs_type=='projected'&& y_name %in% 
     c('lat', 'latitude','Lat', 'Latitude')){
    warning(paste0('Mismatch between nc dimnames and expected format, given selected crs type: ', 
                   crs_type, ', y_name:', y_name)
    )}
  if(crs_type == 'projected'){
  ncatt_put(nc, x_name, "standard_name", "projection_x_coordinate")
  ncatt_put(nc, x_name, "long_name", "Easting")
  ncatt_put(nc, x_name, "axis", "X")
  
  ncatt_put(nc, y_name, "standard_name", "projection_y_coordinate")
  ncatt_put(nc, y_name, "long_name", "Northing")
  ncatt_put(nc, y_name, "axis", "Y")
  }else if(crs_type == 'geographic'){
        ncatt_put(nc, "lon", "standard_name", "longitude")
    ncatt_put(nc, x_name, "long_name", "longitude")
    ncatt_put(nc, x_name, "axis", "X")
    
    ncatt_put(nc, y_name, "standard_name", "latitude")
    ncatt_put(nc, y_name, "long_name", "latitude")
    ncatt_put(nc, y_name, "axis", "Y")
  }
}
