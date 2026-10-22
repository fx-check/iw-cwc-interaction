#' Detrends a DEM with `lm.fit()`
#' @param dem a `RasterLayer`
#' @return the detrended raster as a `RasterLayer`
#' @export
#' @keywords fft2d
detrend_dem <- function(dem){
  coord <- raster::coordinates(dem)
  x <- coord[, 1]
  y <- coord[, 2]
  z <- raster::getValues(dem)
  ind <- which(!is.na(z))
  
  # odr <- pracma::odregress(coord[ind, ], z[ind])
  # z_detrend[ind] <- z[ind] - odr$fitted
  lmr <- stats::lm.fit(coord[ind, ], z[ind])
  z_detrend <- z
  z_detrend[ind] <- lmr$residuals
  
  detrended_dem <- dem
  detrended_dem <- raster::setValues(detrended_dem, z_detrend)
  detrended_dem[is.na(detrended_dem)] <- 0
  return(detrended_dem)
}

#' Detrends a matrix dem with `lm.fit()`
#' @param mat a matrix from `RasterLayer`
#' @param coord the coordinates from the associated raster
#' @return the detrended matrix
#' @export
#' @keywords fft2d
detrend_mat <- function(mat, coords){
  # coord <- raster::coordinates(mat)
  # x <- coords[, 1]
  # y <- coords[, 2]
  # z <- raster::getValues(mat)
  ind <- which(!is.na(mat))
  
  # odr <- pracma::odregress(coord[ind, ], mat[ind])
  # z_detrend[ind] <- mat[ind] - odr$fitted
  lmr <- stats::lm.fit(coords[ind, ], mat[ind])
  z_detrend <- mat
  z_detrend[ind] <- lmr$residuals
  z_detrend[is.na(z_detrend)] <- 0
  # detrended_dem <- mat
  # detrended_dem <- raster::setValues(detrended_dem, z_detrend)
  # detrended_dem[is.na(detrended_dem)] <- 0
  return(z_detrend)
}

#' Pads a matrix with zeros so it has a given dimension
#' 
#' This function is useful to pad a matrix so that the dimension is a factor of 2.
#' 
#' @param x a `matrix`
#' @param n the desired number of rows of `x`
#' @param m the desired number of columns of `x`
#' @return a padded version of `x` with dimension `n,m`
#' @export
#' @keywords fft2d
pad_mat <- function(x, n, m){
  padx <- n - nrow(x)  
  pady <- m - ncol(x) 
  padx <- padx / 2
  pady <- pady / 2
  x <- matlab::padarray(x, c(floor(padx), floor(pady)), direction = "pre")
  x <- matlab::padarray(x, c(ceiling(padx), ceiling(pady)), direction = "post")
  return(x)
}

#' A two-dimensional Hann window
#' 
#' This code is (largely) a port from Matlab to `R` of the code from [2DSpecTools](https://web.mit.edu/perron/www/files/2DSpecTools-v1.1.zip): A 2D spectral analysis toolkit for Matlab by Taylor Perron et al. 
#' 
#' @param M a `matrix`
#' @return a list with two elements the windowed matrix `H` and the sum of the coefficients `Wss`
#' @export
#' @keywords fft2d
Hann2D <- function(M){
  nx <- ncol(M)
  ny <- nrow(M)
  # matrix coordinates of centroid of M
  a <- (nx + 1) / 2 
  b <- (ny + 1) / 2 
  XY <- pracma::meshgrid((1:nx), (1:ny))
  X <- XY$X
  Y <- XY$Y
  theta <- (X == a) * (pi / 2) + (X != a) * atan2((Y - b), (X - a))
  r <- sqrt((Y - b)^2 + (X - a)^2)
  rprime <- sqrt((a^2) * (b^2) * (b^2 * (cos(theta))^2 + a^2 * (sin(theta))^2)^(-1))
  hanncoeff <- (r < rprime) * (0.5 * (1 + cos(pi * r / rprime)))
  H <- M * hanncoeff
  Wss <- sum(sum(hanncoeff^2))
  return(list(H = H, Wss = Wss, hanncoeff = hanncoeff))
}



#' Perfoms a two-dimensional Fourier transform
#' 
#' This code is (largely) a port from Matlab to `R` of the code from [2DSpecTools](https://web.mit.edu/perron/www/files/2DSpecTools-v1.1.zip): A 2D spectral analysis toolkit for Matlab by Taylor Perron et al. 
#' As `dx` and `dy` refers to spacing in map units, it is recommended, if applicable, that the original raster is projected from latitude, longitude to a planar projection.
#' 
#' @param M a `matrix`
#' @param dx the spacing in the x direction in map units (usually meters) 
#' @param dy the spacing in the y direction in map units (usually meters)
#' @param dr the spacing in x & y direction in map units (usually meters), overrides `dx` and `dy` if provided
#' @param Hann `logical`, if `TRUE` performs a Hann windowing, default to `TRUE`
#' @return a list with 8 elements: `spectral_power_matrix`, `radial_frequency_matrix`, `spectral_power_vector`, `1/L`, `radial_frequency_vector`, `wave_vectork_xy`.
#' `kx` and `ky`, wave vector matrices matching the spectral_power_matrix from -2π×dr/L to 2π×dr/L (in rad/unit of `dr`), including K=0 offset from the centre by 1 index position in the positive.
#' The wave vectors are identical when dr == dx == dy, otherwise the input matrix is padded to the larger dimension 
#' The `spectral_power_matrix` is the window-corrected 2d DFT periodogram returned in units of dr^2, typically m^2, following normalisation by the window function
#' The  `spectral_power_vector` is the 1d DFT periodogram along monotonically increasing frequencies, multiplied by two to account for returning only half the spectrum
#' `dfx` and `dfy` is frequency (e.g. 1/(dx*L)) in cycles per unit length of `dx`/`dy`/`dr`, typically meters, multiply by 2π for angular wavenumber
#' `radial_frequency_matrix` and `radial_frequency_vector` are provided in (angular) radial wavenumbers (rad/input unit)
#' @export
#' @keywords fft2d
fft2D <- function(M, dx, dy, dr = NULL, Hann = TRUE){
  nx <- ncol(M)
  ny <- nrow(M)
  if(!is.null(dr)){dx = dy = dr} #dr overrides dx & dy if provided
  stopifnot('dx and dy or dr must be provided and positive'= 
              (all(!is.null(dx),!is.null(dy)) | !is.null(dr)) & all(c(dx,dy,dr)>0))
  if(abs(nx-ny)/max(nx,ny)>0.2) warning('The supplied matrix deviates from square dimensions by >20%, 
                                        proceeding but suggesting caution')
  if (Hann){
    l <- Hann2D(M)
    M <- l$H
    Wss <- l$Wss
  } else {
    Wss <- sum(sum(pracma::ones(ny, nx)))
  }
  Lx <- 2^(ceiling(log(max(c(nx, ny)))/log(2))) # mandates even Lx, important below
  Ly <- Lx
  dfx <-  1 / (dx * Lx) # frequency in cycles per unit(dx)
  dfy <-  1 / (dy * Ly)
  dk <- 2 * pi / (dx * Lx) # wave vector spacing
  M <- pad_mat(M, Lx, Ly) #DC at index Lx/2 +1, Ly/2 +1
  # replicate in wave vectors below
  k_x <- (-Lx/2):(Lx/2 - 1) * (2*pi*dfx) # multiplication of index by angular wavenumber 
  k_y <- (-Ly/2):(Ly/2 - 1) * (2*pi*dfy) # produces wavenumber vector for every fft cell
  #the resulting wave vectors are in units of radians per unit of xy, typcially meters
  #matrices from k_x and k_y:
  KX <- matrix(rep(k_x, each = Ly), nrow = Ly, ncol = Lx)
  KY <- matrix(rep(k_y, times = Lx), nrow = Ly, ncol = Lx)
  
  M <- stats::fft(M)
  M <-  mrbsizeR::fftshift(M)
  M[(Ly/2 + 1), (Lx/2 + 1)] <- 0
  #Calculate the DFT periodogram, which has units of amplitude^2(m^2 for 
  # topography). Dividing by (Wss*Lx * Ly) corrects for the reduction of
  # amplitude by the windowing function, with Wss == Lx*Ly if no windowing used
  M <- M * Conj(M) / (Lx * Ly * Wss)
  spectral_power_matrix <- M 
  xc <- Lx / 2 + 1 # matrix indices of zero frequency
  yc <- Ly / 2 + 1
  XY <- pracma::meshgrid((1:Lx), (1:Ly)) # matrices of column and row indices
  cols <- XY$X
  rows <- XY$Y
  
  radial_frequency_matrix <- sqrt((dfy * (rows - yc))^2 + (dfx*(cols - xc)) ^2) * 2 * pi
  # Creating sorted, non-redundant vectors of frequency and power 
  M <- M[ , (1:(Lx/2+1))]
  radial_frequency_vector <- radial_frequency_matrix[, (1:(Lx/2+1))]
  #' The data is symmetric around zero. Set the frequencies to negative values,  
  #' they will be clipped out further down
  radial_frequency_vector[((yc+1):Ly), xc] <- -1  
  # concatenate column vectors of frequency and PSD and sort by frequency:
  radial_frequency_vector <- pracma::sortrows(cbind(c(radial_frequency_vector), 
                                                    c(Re(M))), 1)
  #' Take only positive frequencies. This gets rid of DC (zero freq) as well as 
  #' the redundant frequencies we set to -1 above:
  radial_frequency_vector <- radial_frequency_vector[radial_frequency_vector[ , 1] > 0 , ]
  # Separate into power and frequency vectors and assign to output arguments; 
  # the factor of 2 corrects for the fact that we have taken only half of the 2D 
  # spectrum_ sum(spectral_power_vectorec) should now equal sum(Pmat(:))_
  spectral_power_vector <-  2 * radial_frequency_vector[ , 2]  
  radial_frequency_vector <- radial_frequency_vector[, 1]
  l <- list(
    spectral_power_matrix = spectral_power_matrix,
    radial_frequency_matrix = radial_frequency_matrix,
    spectral_power_vector = spectral_power_vector,
    radial_frequency_vector = radial_frequency_vector,
    kx = KX,
    ky = KY,
    dfx = dfx,
    dfy = dfy)
  return(l)
}

#' Tranforms kx, ky matrices of wave vectors to be aligned with flow u,v


#' @param `kx`,`ky` wave vector matrices returned by fft2d, or `fft2d.out` ouptput object that overrides `kx` and `ky`
#' @param `u`,`v` flow vector components in the eastward and northward direction respectively
#' @return A list with flow adjusted wave vector matrices `k_along` and `k_across`, as well as a named vector `flow` with flow direction as CCW from east in radians (used for wave vectors), flow magnitude and flow direction CW from north (in degrees) for geographic standards
#'  
align_k <- function(kx, ky, u, v, fft2D.out = NULL){
  if(!is.null(fft2D.out)){
    kx <- fft2D.out$kx
    ky <- fft2D.out$ky
  }
  stopifnot('fft2d.out invalid or kx and ky need to be supplied as matrices of the same dimension'=
              all(is.matrix(kx), is.matrix(ky),
                  all(dim(kx) == dim(ky))))
  stopifnot('u and v need to be supplied as scalar values'= 
              all(length(u)==1L,
                  length(v)==1L))
  theta <- atan2(v, u) #flow direction from east, CCW
  k_along <- kx * cos(theta) + ky * sin(theta)
  k_across <- ky * cos(theta) - kx * sin(theta)
  if(is.null(fft2D.out)) {k_mag <- sqrt(kx^2+ky^2)} else {k_mag <- fft2D.out$radial_frequency_matrix} 
  k_mag_flow <-  sqrt(k_along^2 + k_across^2)
  phi <- atan2(u, v)*180/pi #flow direction from north, CW in degrees
  flow_mag <- sqrt(u^2 + v^2)
  if (!all(dplyr::near(k_mag- k_mag_flow,0))) warning('Mismatch in wave vector magnitude calculations')
  flow <- c('theta_CCW_x_rad' = theta,
            'phi_CW_N_degree' = phi,
            'magnitude' = flow_mag)
  l <- list(k_along = k_along,
            k_across = k_across,
            flow = flow)
  return(l)
}

