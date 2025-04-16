library(roxygen2)
actd_spline <- function(ctd, x_var = 'pressure', y_var = 'sigma0', min_d, bin_size, envir = environment()){
  #' generates adaptive splines based on inflection points of the y_var, typicall sigma for N2 calculation 
  #' if x_var is 'pressure' and y_var is 'sigma' density input (note columns need to be named accordingly), smoothed N2 profiles are generated
  #' @param ctd a tibble or data frame of ctd, typically with sigma (pot density) and pressure data
  #' @param x independent data for splines, typically pressure (or depth). Supply variable name (as in column name) as string
  #' @param y dependent data variable for splines, typically sigma0 (though different ref pressure may be chosen). Supply variable name (as in column name) as string
  #' @param min_d minimum x value for spline
  #' @param bin_size bin size (search distance) for second order differencing -> finding inflection knots
  #' @param envir choose the environment in which the splines will be stored
  library(tidyverse)
  x_dat <- get(x_var,pos = ctd)
  y_dat <- get(y_var,pos = ctd)
  cast <- list() # list as container for all data
  cast$ctd <- ctd
  p_int <- unique(diff(x_dat))
  if(length(p_int)>1){return('The CTD supplied has irregular pressure sampling intervals, please filter or bin it prior to running this function')}
  l <- bin_size/p_int # 50 m lags, taking into account the  pressure interval
  if((bin_size %% p_int) !=0){return(sprintf("The bin size is not divisible by the pressure interval %.1f, please supply an appropriate and divisible bin size", p_int))}
  
  max_d <- max(x_dat, na.rm = T)
  xnew = seq(min_d, max_d, by = p_int)
  
  # 2nd order of difference, here 2nd order - i.e. looking for inflections
  ds = diff(y_dat, lag = l, differences = 2) #2nd order diffs over bin_size
  
  neigh_bool <- dplyr::near(ds,0, tol = 0.0005) # potential knots are dplyr::near zero, i.e. inflection points
  pt_knots <- x_dat[neigh_bool] # select knots based on boolean neighbourhood vector
  pt_knots <- pt_knots[pt_knots>min_d] # remove knots below minimum spline depth 
  
  ### Below, we group knots to remove the ones that are too close together!
  abc = 1
  trsh = 6 # minimum distance for knots to be kept for B-spline
  # p.int doesn't matter here because we are now talking pressures, not lags
  diff_knot_groups <- c(LETTERS[abc]) # assign group 'A' to the first knot
  diffs <- diff(pt_knots, lag = 1) # assess distance between knots
  for (i in 2:length(pt_knots)) {
    if (diffs[i-1]<=trsh){ #if distance between knots is <6 m, keep in group
      diff_knot_groups <- append(diff_knot_groups, diff_knot_groups[i-1])
    } else if (abc<length(LETTERS)) { # if distance greater than trsh, put in new group
      # with single letter from alphabet unless we get to Y (end of alphabet)
      abc = abc+1
      diff_knot_groups <- append(diff_knot_groups, LETTERS[abc])
    }
    else  { # once alphabet has been used, below code extends this to two letters
      # providing up to 26^2 == 676 groups
      abc = abc+1
      ix = floor(abc/length(LETTERS))
      jx = abc%%length(LETTERS)+1 # plus one to avoid 0
      diff_knot_groups <- append(diff_knot_groups, 
                                 paste0(LETTERS[ix],LETTERS[jx]))
    }
  }
  diff_knots <- sort(aggregate(pt_knots, by = list(diff_knot_groups), mean)$x)
  diff_knots_b <- numeric()
  j=0
  # add more knots when distance between inflection points exceeds multiples of 
  # the bin size:
  for (i in 1:(length(diff_knots)-1)){
    diff_knots_b <- append(diff_knots_b, diff_knots[i], after = j)#after not at
    j = j+1
    if (diff_knots[i]+1.9*bin_size < diff_knots[i+1]){
      a <- diff_knots[i]
      b <- diff_knots[i+1]
      #regular spaced new knots including a & b
      additional_knots <- seq(a, b, by = (b-a)/ceiling((b-a)/bin_size))
      #exclude a & b: 
      additional_knots <- additional_knots[-c(1,length(additional_knots))]
      #append to diff_knots_b, rinse & repeat, going onto the next knot distanc check
      diff_knots_b <- append(diff_knots_b,additional_knots, after = j)
    }
  }
    
    # Generate knots at start if no inflection points found
  pre_kts <- seq(min(xnew),
                 diff_knots_b[1],
                 by = bin_size)[-1]
  diff_knots_b <- append(diff_knots_b, 
                         values = pre_kts[-length(pre_kts)],
                         after = 0)
  post_kts <- c() # knots after last diff knot in case no inflection points found in th last bin(s)
  if (diff_knots_b[length(diff_knots_b)]<(x[length(x)]-bin_size)){
    post_kts <- seq(diff_knots_b[length(diff_knots_b)],max(xnew), by = l)
    
  }
  if (!is.null(post_kts)){
    if (post_kts[length(post_kts)]==xnew[length(xnew)]) {
      post_kts <- post_kts[-length(post_kts)]} # in the rare event that this 
    # sequence ends in a terminal kt it is removed here, otherwise we get 
    # rank-deficient predictions
  } #The outer if checks that there are post_kts at all, avoid length zero error
  cast$dkts <- append(diff_knots_b,
                     values = post_kts[-1])
  # BM <- splines::bs(x = xnew, knots = cast$dkts) # fishy, error in b.sp
  BM <- splines::bs(x = x_dat, knots = cast$dkts)
  b.sp <- lm(y_dat~., data = as.data.frame(BM))
  cast$b.sp <- b.sp$call # for tracebility, whole model is messy to store in list
  newBM = predict(BM, xnew)
  
  cast$pred.bspline = as.data.frame(cbind(x_var= xnew, 
                                         y_var = as.vector(predict(b.sp, 
                                                                      newdata=as.data.frame(newBM),
                                                                      rankdeficient = 'NAwarn')
                                         )
  )
  )
  if((x_var=='pressure' | x_var=='Pressure' | x_var=='p') &
     (startsWith(y_var, 'sigma') | startsWith(y_var, 'Sigma'))
  ){
    cast$N2 <- as_tibble(cbind('pressure'=xnew, 
                               'N2' = oce::swN2(pressure = xnew, 
                                                sigmaTheta = get('y_var',cast$pred.bspline))
    ))  
  }
  
  return(cast)
}
