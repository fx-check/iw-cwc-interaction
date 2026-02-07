
# Here, we construct a function to compute a bottom weights array for 
# depth-weighted bottom layer averaging based on true depth and HYCOM 
# grid depth levels. The output is a 3D array that may be transformed to a stars
# object.
# Function defininition ---------------------------------------------------


build_n2_bottom_weights <- function(H_true, max_n2_depth, depth_levels, force_last_bin =TRUE,
                                    Lfun = function(H) {pmin(pmax(0.25 * H, 300), 1200)}) {
  #' H_true: matrix of true depth
  #' H_max_level: matrix of max depth level of HYCOM grid, same dims as H_true
  #' depth_levels: vector of depth levels from N2 data with length K. Typically these are 
  #' mid depths between observed/modelled layers, by nature of finite differencing 
  #' Lfun: function to compute bottom averaging window thickness L based on depth H
  #' Returns: 3D array of bottom weights, dim (nx, ny, nz=K)
  #' This specific implementation is targeted at weighing N2 data
  #' 
  if(inherits(depth_levels, 'units')) depth_levels  <- units::drop_units(depth_levels)
  K <- length(depth_levels)
  n2_edges <-  c(0, depth_levels[-K] +diff(depth_levels)/2,
                 depth_levels[K] +diff(depth_levels)[K-1]/2) 
  #side idea: RUN A SPLINE AND DEPTH INTEGRATE BASED ON THE SPLINE FUNCTION or use S-G filter!
  #PCHIP / monotone cubic (R: splinefun(method="monoH.FC", Barker & McDougall:
  # https://journals.ametsoc.org/view/journals/atot/37/4/jtech-d-19-0211.1.xml
  
  nx <- dim(H_true)[1]; ny <- dim(H_true)[2]
  stopifnot('H_true and max_n2_depth must have same dimensions' = all(dim(H_true)==dim(max_n2_depth)))
  w <- array(0, dim = c(nx, ny, K))
  n2_top <- n2_edges[-length(n2_edges)] #exclude deepest edge from N2 tops
  n2_bot <- n2_edges[-1]    #exclude shallowest edge from N2 bottoms
  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      H <- H_true[i, j] # true depth
      if (is.na(H) || H <= 0) next  # land or invalid
      z_max <- max_n2_depth[i,j] # max N2 mid-depth at this location
      if (!is.finite(z_max)) next
      #identify depth bin in which to locate z_max:
      z_idx <- findInterval(z_max, n2_edges, left.open = TRUE) # 1..nz
      z_idx <- max(1L, min(z_idx, K)) #safe guard
      
      L <- Lfun(H)
      if (L> H){L<-H} #reduces L for water 'shallower' than L
      z_0 <- max(0, H - L)  # top (in water column) of averaging window
      z_1 <- min(H, n2_edges[z_idx + 1L]) # max depth to which we want to depth-average, which 
      # is the shallower of the true depth and the deepest valid N2 cell *edge*
      #put all weights in the last valid N2 bin if z_0 deeper than z_1:
      if (z_1<=z_0) {
        if (force_last_bin) {w[i,j,] <- 0; w[i,j,z_idx] <- 1}
        # else leave as zeros 
        next
      }
      # Fast overlap with each cell: the base extremes recycles z_1 and z_0, 
      # comparing them to each element of n2_bot and n2_top respectively,
      # this isolates the n2 'cells' [top, bottom] that overlap with [z_0, z_1] window 
      
      overlap_n2 <- pmax(0, #zero unless the below subtraction is positive
                         pmin(n2_bot, z_1) - #n2 bottom down to z_1 max averaging depth
                           pmax(n2_top, z_0)) #z0 down to top of window, then each n2_top
      L_eff <- sum(overlap_n2) #like z_1 - z_0: effective overlap in case H>zmax, ==L where H<z_max
      # scale by effective L to get weights [0,1]:
      w[i, j, ]  <- if (L_eff > 0) overlap_n2 / L_eff else rep(0, length(overlap_n2))
      
    }
  }
  return(w) 
}
