utm18n <- "EPSG:32618"  # UTM zone 18N

#--- Helper: grid convergence gamma (radians) at a lon/lat point -----------
# Compute local angle between grid north (UTM y-axis) and true/geodetic north.
# Method: step ~100 m along true north & east on the sphere, project both to UTM,
# then measure their directions in projected (x,y).
# refer to IOGP for guidance on conventions used
 # https://www.iogp.org/bookstore/wp-content/uploads/sites/2/woocommerce_uploads/2017/01/373-21.pdf
gamma_utm <- function(lon, lat, utm_crs = utm18n, step_m = 100,#new below:
                      convention = c('Gauss-Bomford',"Survey"), flip_SH = F
){
  #' @param convention  "Gauss-Bomford": subtracts grid convergence from azimuth to
  #' obtain bearing, convergence is positive when True North lies west of Grid North
  #' "Survey": adds grid convergence from azimuth to obtain bearing, convergence 
  #' is positive when True North lies east of Grid North
  convention <- match.arg(convention)
  central_meridian <- mean(crs(utm18n, describe=T)$extent[[1]][1:2])
  if(abs(lon-central_meridian)>6){
    warning(cat('Distance from central meridian at', central_meridian, 
                'to chosen longitude at', lon, 'exceeds 6°', sep = " "))}
  p0_ll <- c(lon, lat)
  pN_ll <- destPoint(p0_ll, b = 0,  d = step_m)  # towards true north
  # project to UTM
  p0_xy <- as.numeric(project(matrix(p0_ll, nrow = 1), "EPSG:4326", utm_crs))
  pN_xy <- as.numeric(project(matrix(pN_ll, nrow = 1), "EPSG:4326", utm_crs))
  vN    <- pN_xy - p0_xy #added: variation in xy through UTM projection
  # the above would default to 0,step_m on central meridian

  
  gamma_raw <- atan2(-vN[1], vN[2])  # radians, defaults to Gauss-Bomford on NH
  if(flip_SH & lat<0){gamma_raw <- -gamma_raw} # flips sign for GB on SH
  gamma <- switch(convention,
                  'Gauss-Bomford' =  gamma_raw,
                  'Survey' = -gamma_raw) #flip signs for survey convention
  names(gamma)<- paste(convention, 'convention')
  return(gamma)
}