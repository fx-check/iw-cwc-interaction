# Load required packages
required_packages <- c('raster','oce', 'stars','tidyverse','cubelyr')
for(i in required_packages) {library(i, character.only = T)}


# Helper function definitions ---------------------------------------------

# The following functions are defined to analyse TS data from HYCOM modelling output


# Check and initialize parameters -----------------------------------------
check_init<- function(x_star, max_depth = NULL){
  # Check if all required parameters are present
  required_vars <- c("salinity", "water_temp")
  missing_vars <- setdiff(required_vars, names(x_star))
  if (length(missing_vars) > 0) {
    stop(paste("Missing required parameters:", paste(missing_vars, collapse = ", ")))
  }
  required_dims <- c("x", "y", "depth")
  missing_dims <- setdiff(required_dims, names(st_dimensions(x_star)))
  if (length(missing_dims)>0){
    stop(paste("Missing required dimensions:", paste(missing_dims, collapse = ", ")))
  }
  # Grab depth levels and trim to max_depth if specified
  depth_levels <- st_get_dimension_values(x_star, 'depth')
  if(inherits(depth_levels, "units") && 
     !identical(units::deparse_unit(depth_levels), "m")){
    warning('Depth units may not be in meters, proceeding under the assumption that depth levels are in meters, positive with increasing depth. Please ensure data fits this expected format.')
  }
  depth_levels <- units::drop_units(depth_levels)
  if (any(diff(depth_levels) <= 0)) {
    stop("Depth levels are not strictly increasing. Please fix ordering before initialization.")
  }
  if (!is.null(max_depth)) {
    # extract deepest level shallower than specified max detph
    idx <- which(depth_levels <= max_depth)
    if (length(idx) == 0L) stop("max_depth is shallower than the shallowest level; nothing to keep.")
    idx_max <- max(idx)
    # slice data to eliminate levels deeper than max depth:
    x_star <- slice(x_star, 'depth', index = 1:idx_max)
    depth_levels <- depth_levels[1:idx_max]
  }
  # Include depth dimension as attribute, this is useful to retain,
  # as we will later convert the vertical dimension to pressure
  x_star$depth <- st_dim_to_attr(x_star, which = which(names(st_dimensions(x_star))=='depth'))
  # mask depth levels where both salinity and temperature are NA, either land or seabed
  combined_mask <- is.na(x_star$salinity) & is.na(x_star$water_temp)
  x_star$depth[combined_mask] <- NA
  # Calculate mean latitude for pressure calculation
  y_vals <- st_get_dimension_values(x_star, "y")
  mean_lat <- mean(y_vals)
  # Convert longitude for later calcs to stars attribute, [0-360] for TEOS:
  
  x_star$longitude <- st_dim_to_attr(x_star, which = which(names(st_dimensions(x_star))=='x')) %% 360
  # Convert latitude for later cals to stars attribute:
  lats <- st_get_dimension_values(x_star, 'y')
  x_star$latitude  <- st_dim_to_attr(x_star, which = which(names(st_dimensions(x_star))=='y'))
  # Calculate pressure levels
  p_levels <- gsw::gsw_p_from_z(-depth_levels, lat = mean_lat)  # in dbar
  # if grids with a larger latitude span were considered, this should be modfified to account for
  # pressure dependence on latitude. Here, this only introduces a 0.025% difference, 
  # which is minute in our depth range
  
  # Rename the vertical dimension to "pressure" and attach p_levels as its coordinate
  x_star <- st_set_dimensions(x_star, "depth", names = "pressure")
  x_star <- st_set_dimensions(x_star, "pressure", values = p_levels)
  
  return(x_star)
}

#cleaned up version, this still has a snag which needs fixing
check.init <- function(x_star, max_depth = NULL, verbose = TRUE) {
  # --- 0) Required fields and dims -------------------------------------------
  required_vars <- c("salinity", "water_temp")
  missing_vars <- setdiff(required_vars, names(x_star))
  if (length(missing_vars) > 0) {
    stop(paste("Missing required parameters:", paste(missing_vars, collapse = ", ")))
  }
  required_dims <- c("x", "y", "depth")
  missing_dims <- setdiff(required_dims, names(st_dimensions(x_star)))
  if (length(missing_dims) > 0) {
    stop(paste("Missing required dimensions:", paste(missing_dims, collapse = ", ")))
  }
  
  # --- 1) Depth vector and checks --------------------------------------------
  depth_levels <- st_get_dimension_values(x_star, "depth")
  
  # Units & conversion
  if (inherits(depth_levels, "units")) {
    if (!identical(units::deparse_unit(depth_levels), "m") && verbose) {
      warning("Depth units are not 'm'; proceeding assuming meters, positive downward.")
    }
    depth_levels <- units::drop_units(depth_levels)
  } else {
    # No units present, surpress this message with verbose=F if undesired
    if(verbose){
      message("Depth levels have no units; assuming meters, positive downward.")
    }
  }
  # Strict monotonicity
  if (any(diff(depth_levels) <= 0)) {
    stop("Depth levels are not strictly increasing. Please fix ordering before initialization.")
  }
  
  # --- 2) Optional: trim to max_depth ----------------------------------------
  if (!is.null(max_depth)) {
    idx <- which(depth_levels <= max_depth)
    if (length(idx) == 0L) {
      stop("max_depth is shallower than the shallowest level; nothing to keep.")
    }
    idx_max <- max(idx)
    x_star <- slice(x_star, "depth", index = 1:idx_max)
    depth_levels <- depth_levels[1:idx_max]
  }
  
  # --- 3) Broadcast depth attribute and mask where both SP & T are NA --------
  x_star$depth <- st_dim_to_attr(x_star, "depth")
  combined_mask <- is.na(x_star$salinity) & is.na(x_star$water_temp)
  x_star$depth[combined_mask] <- NA_real_
  
  # --- 4) Mean latitude and lon/lat attributes -------------------------------
  y_vals <- st_get_dimension_values(x_star, "y")
  mean_lat <- mean(y_vals)
  #usually longitude is already in [0, 360) for HYCOM, but just in case:
  x_star$longitude <- st_dim_to_attr(x_star, "x") %% 360  # [0, 360) for TEOS-10
  x_star$latitude  <- st_dim_to_attr(x_star, "y")
  
  # --- 5) Compute sea pressure from depth (dbar) -----------------------------
  p_levels <- gsw::gsw_p_from_z(z = -depth_levels, lat = mean_lat)
  
  # --- 6) Rename vertical dim to 'pressure' and attach coordinates -----------
  x_star <- st_set_dimensions(x_star, "depth", names = "pressure")
  x_star <- st_set_dimensions(x_star, "pressure", values = p_levels)
  
  return(x_star)
}



 
# Absolute salinity calculation -------------------------------------------
ts.sa <- function(x_star){
  # --- 1) Checks: ------------------------------------------------------
  required_vars <- c("salinity", "water_temp", "longitude", "latitude")
  missing_vars <- setdiff(required_vars, names(x_star))
  if (length(missing_vars) > 0) {
    stop(paste("Missing required parameters:", paste(missing_vars, collapse = ", "),
               'Have you rund initialization function first?'))
  }
  required_dims <- c("pressure", "time")
  missing_dims <- setdiff(required_dims, names(st_dimensions(x_star)))
  if (length(missing_dims) > 0) {
    stop(paste("Missing required dimensions:", paste(missing_dims, collapse = ", "),
               'Have you rund initialization function first?'))
  }
  # --- 2) Pre-allocate absolute salinity array and get dims ---------------------
  SA_arr <- array(NA_real_, dim = dim(x_star))  # x × y × pressure * time
  K <- st_dimensions(x_star)$pressure$to - st_dimensions(x_star)$pressure$from +1L
  trange <- st_dimensions(x_star)$time$to - st_dimensions(x_star)$time$from +1L
  stopifnot('Pressure or time dimension is unaligned with position 3&4 expected here'=
              dim(x_star)[3]==K & dim(x_star)[4]== trange)
  # --- 2) Calculate absolute salinity and assign to array ---------------------
  # calculate absolute salinity for all timesteps, depth levels and xy cells:
  for (i in seq_len(trange)){
    for (k in seq_len(K)) {
      SP_k <- x_star[['salinity']][,,k,i]  # PSU for depth level slice
      SA_arr[,,k,i] <- gsw::gsw_SA_from_SP(
        SP        = SP_k,
        p         = p_levels[k],   # scalar for this slice (in dbar)
        longitude = x_star$longitude[,,k,i],
        latitude  = x_star$latitude[,,k,i]
      )
    }
  }
  # Store back into stars
  SA_stars <- st_as_stars(SA_arr, dimensions = st_dimensions(x_star)) %>%
    setNames("SA")
  x_star$SA <- SA_stars$SA
  return(x_star) 
}
# Conservative temp computation for stars data ----------------------------

ts.CT <- function(x_star){
  # --- 1) Checks: ------------------------------------------------------
  required_vars <- c("SA", "water_temp")
  missing_vars <- setdiff(required_vars, names(x_star))
  if (length(missing_vars) > 0) {
    stop(paste("Missing required parameters:", paste(missing_vars, collapse = ", "),
               'Have you run initialization and ts.sa() functions first?'))
  }
  required_dims <- c("pressure", "time")
  missing_dims <- setdiff(required_dims, names(st_dimensions(x_star)))
  if (length(missing_dims) > 0) {
    stop(paste("Missing required dimensions:", paste(missing_dims, collapse = ", "),
               'Have you rund initialization function first?'))
  }
  # --- 2) Pre-allocate absolute salinity array and get dims ---------------------
  
  # pre-allocate CT array and get pressure levels:
  CT_array <- array(NA_real_, dim = dim(x_star))  # x × y × pressure * time
  p_levels <- st_get_dimension_values(test2, 'pressure')
  trange <- st_dimensions(x_star)$time$to - st_dimensions(x_star)$time$from +1L
  # add data to array:
  for (i in seq_len(trange)){
    for(p in seq_along(p_levels)){
      CT_array[,,p,i] <- gsw_CT_from_t(t = x_star[['water_temp']][,,p,i],
                                     SA = x_star[['SA']][,,p,i],
                                     p = p_levels[p])
    }  
  }
  CT <- st_as_stars(CT_array, dimensions = st_dimensions(x_star)) %>%
    setNames("CT")
    # units(CT$CT) <- as_units("degree_C")
  #add this CT array as an attribute to the TS array, check how it deals with time steps:
  x_star$CT <- CT
  return(x_star)
}
  
# N2 computation for stars data ----------------------------------------
# per profile N2 calculation helper function to hanlde NAs
compute.N2.profile <- function(SA, CT, p, lat) {
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
# bottom N2 mask
depth_levels <- -round(gsw_z_from_p(p_levels, 31.5),0)
depth_mask <- array(data= FALSE, dim = dim(test2))
for (n in seq_len(nx)) {
  for (m in seq_len(ny)) {
    profile <- test2$depth[n,m,,]
    profile_mask <- !is.na(profile)
    max_depth_idx <- which.max(profile) #insusceptible to NAs
    max_depth <- max(profile, na.rm = T)
    if(max_depth>2500) {# only use one 1km N2 value for bottom N2 below 2500m
      profile_mask[-max_depth_idx] <- FALSE} # where spacing depth levels are 1km 
    else if(max_depth> 1250){# use 2 ~500m spaced levels for 1500-2500m
      # if seabed depth is say 1600m, this will also use 1250 level
      # reducing averaging window as water depth gets shallower
      profile_mask[-c((max_depth_idx-1):max_depth_idx)] <- FALSE
    }
    else if(max_depth> 350){# use 3 100-250m spaced levels for 400-1250m
      # reducing averaging window further as water depth gets shallower
      profile_mask[-c((max_depth_idx-2):max_depth_idx)] <- FALSE
    }
    #otherwise, use all values to surface
  }
}


ts.N2 <- function(x_star) {
  # --- 1) Checks and grab required (dim) values ------------------------------
  required_vars <- c("SA", "CT", 'latitude')
  missing_vars <- setdiff(required_vars, names(x_star))
  if (length(missing_vars) > 0) {
    stop(paste("Missing required parameters:", paste(missing_vars, collapse = ", "),
               'Have you run initialization, ts.sa() and ts.ct() functions first?'))
  }
  required_dims <- c("pressure", "time")
  missing_dims <- setdiff(required_dims, names(st_dimensions(x_star)))
  if (length(missing_dims) > 0) {
    stop(paste("Missing required dimensions:", paste(missing_dims, collapse = ", "),
               'Have you rund initialization function first?'))
  }
  K <- st_dimensions(x_star)$pressure$to - st_dimensions(x_star)$pressure$from +1L
  trange <- st_dimensions(x_star)$time$to - st_dimensions(x_star)$time$from +1L
  nx <- st_dimensions(x_star)$x$to - st_dimensions(x_star)$x$from +1L
  ny <- st_dimensions(x_star)$y$to - st_dimensions(x_star)$pressure$from +1L
  
  # --- 2) Compute N2 array ------------------------
  for (i in seq_len(trange) ) {
    for (n in seq_len(nx)) {
      for (m in seq_len(ny)) {
        lat <- x_star$latitude[1,m,1,1] #1s are placeholders, won't vary along those dims
        SA_col <- x_star[['SA']][n, m, ,i]
        CT_col <- xstar[['CT']][n, m, ,i]
        
        N2_array[n, m, ] <- compute.N2.profile(SA_col, CT_col, p_levels, lat)
      }
    }
  }
  # Wrap back into stars and set the pressure-mid coordinate
  dims_N2  <- dims_xyk
  dims_N2$pressure$values <- p_mid
  dims_N2$pressure$to     <- length(p_mid)
  
  N2_stars <- st_as_stars(N2_array, dimensions = dims_N2)
  names(N2_stars) <- "N2"   # s^-2 (you can add units if you use 'units' package)
  return(N2_stars)
}

# dims --------------------------------------------------------------------

ts_xyk <- st_redimension(ts, new_dims = st_dimensions(ts)[1:3])
dims_xyk <- st_dimensions(ts_xyk)
nx <- dims_xyk[1]
ny <- dims_xyk[2]
K  <- dims_xyk[3]
stopifnot(K == length(p_levels))

