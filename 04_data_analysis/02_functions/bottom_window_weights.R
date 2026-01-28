
# Here, we construct a function to compute a bottom weights array for 
# depth-weighted bottom layer averaging based on true depth and HYCOM 
# grid depth levels. The output is a 3D array that may be transformed to a stars
# object.
# Function defininition ---------------------------------------------------


build_bottom_weights <- function(H_true, H_max_level, depth_levels, 
                                 Lfun = function(H) {pmin(pmax(0.25 * H, 300), 1200)}) {
  #' H_true: raster of true depth
  #' H_max_level: raster of max depth level of HYCOM grid, same dims as H_true
  #' depth_levels: vector of depth level mid depths, length nz
  #' 
  dmid <- diff(depth_levels) #computes inter depth level differences
  # here we assume 0 top edge (surface) and bottom edge is last mid + half last diff:
  depth_edges <- c(0, head(depth_levels, -1) + dmid/2, tail(depth_levels, 1) + tail(dmid, 1)/2)
  dz <- diff(depth_edges) #all diffs with length(dz) equal to depth_levels based on above assumptions
  nx <- dim(H_true)[1]; ny <- dim(H_true)[2]; nz <- length(depth_levels)
  w <- array(0, dim = c(nx, ny, nz))
  mid <- (depth_edges[-1] + depth_edges[-length(depth_edges)])/2
  top <- depth_edges[-(nz+1)]
  bot <- depth_edges[-1]
  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      H <- H_true[i, j] # true depth
      if (is.na(H) || H <= 0) next  # land or invalid
      z_max <- H_max_level[i,j] # max HYCOM depth at this location
      z_idx <- which.min(abs(z_max - mid))
      if (!is.finite(z_max) || !is.finite(mid[z_idx])) next
      
      z_max_edge <- depth_edges[z_idx+1] #associated bottom edge
      L <- Lfun(H)
      if (L> H){L<-H} #reduces L for water 'shallower' than L
      z0 <- max(0, H - L)  # bottom of averaging window
      
      # Fast overlap with each cell [edge[k], edge[k+1]]

      # scale by L to get weights [0,1]:
      overlap <- pmax(0, pmin(bot, H) - pmax(top, z0))/L
      # if the true max depth is deeper than the deepest valid HYCOM cell,
      # add the 'too deep' weigthing to the last valid cell and zero invalid cells:
      overlap[z_idx] <- overlap[z_idx] + sum(overlap[-(1:z_idx)])
      overlap[-(1:z_idx)]<- 0
      
      
      # scaling by L provides weights [0,1] for depth weighted average
      w[i, j, ] <- overlap
    }
  }
 return(w) 
}

# Execution / test --------------------------------------------------------

bottom_weights <- build_bottom_weights(h_nrl_res, max_hycom_depth, depth_edges, dz)
check[check==0] <- NA
View(check)                                              
all(check, na.rm = T)
# TRUE