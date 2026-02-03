# Load DEM tiles from grid and perform FFT with Hann filter
# Write the FFT output including wavenumber vector matrix (for subsequent masking)
# to a dedicated environment with nx ny index. This allows relating the data with
# .nc u&v values for rotation, as well as bottom N2 for energy calculation

source('04_data_analysis/02_functions/fft2d_mod.R')
library(raster)
fft_env <- new.env()
tile_dir <- "C:\\Users\\fbutschek\\Desktop\\hycom\\Bathy_tiles"
tile_names <- list.files(tile_dir, pattern = '.tif$', full.names = F)

# Test run ----------------------------------------------------------------
# tiles have more columns (nc) than rows (nr) due to unequal nc spacing and
# longitude effects
(tile <- raster(paste(tile_dir,tile_names[1200],sep = '/')))
ndif <- diff(dim(tile))[1] #difference between tile dimensions
ndiff_split <- dim(tile)[1] - dim(tile)[2]/2 
#splitting makes the second dimension shorter (nc<nr), let's check wich discrepancy is greater
if( ndif/min(dim(tile)[1:2]) <= ndiff_split/ max(dim(tile)[1:2])/2){ 
  #just calculate one FFT if unsplit diff smaller or equal
  dem_mat <- as.matrix(tile)
  f.out <- fft2D(dem_mat, dr = res(tile)[1], Hann = T)
  nx <- stringr::str_split_i(tile_names[1], pattern = '_', i = 3)
  ny <- stringr::str_split_i(
    stringr::str_split_i(tile_names[1], pattern = '_', i = 5), pattern = '\\.',1)
  fft_name <- paste('nx',nx,'ny',ny,'01', sep = '_')
  assign(fft_name, f.out, envir = fft_env)
} else{ # split in two DEMs and generate 2 FFTs
  n_split <- floor(dim(tile)[2]/2)
  dem_mat <- as.matrix(tile)[,1:n_split]
  f.out <- fft2D(dem_mat, dr = res(tile)[1], Hann = T)
  nx <- stringr::str_split_i(tile_names[1], pattern = '_', i = 3)
  ny <- stringr::str_split_i(
    stringr::str_split_i(tile_names[1], pattern = '_', i = 5), pattern = '\\.',1)
  fft_name <- paste('nx',nx,'ny',ny,'01', sep = '_')
  assign(fft_name, f.out, envir = fft_env)
  #repeat for second split
  dem_mat <- as.matrix(tile)[,(n_split+1):dim(tile)[2]]
  f.out <- fft2D(dem_mat, dr = res(tile)[1], Hann = T)
  fft_name <- paste('nx',nx,'ny',ny,'02', sep = '_')
  assign(fft_name, f.out, envir = fft_env)
}

# Loop across all tiles ---------------------------------------------------
no_split = no_unsplit = 0 # record keeping of unsplit vs split matrices
out.dir <- r'(C:\Users\fbutschek\Desktop\hycom\fft_output)'
for (j in 1:length(tile_names)) {
  tile <- raster(paste(tile_dir,tile_names[j],sep = '/'))
  ndif <- diff(dim(tile))[1] #difference between tile dimensions
  ndiff_split <- dim(tile)[1] - dim(tile)[2]/2 
  
  if( ndif/min(dim(tile)[1:2]) <= ndiff_split/ max(dim(tile)[1:2])/2){ 
    #just calculate one FFT if unsplit diff smaller or equal
    no_unsplit <- no_unsplit+1
    dem_mat <- as.matrix(tile)
    f.out <- fft2D(dem_mat, dr = res(tile)[1], Hann = T)
    nx <- stringr::str_split_i(tile_names[j], pattern = '_', i = 3)
    ny <- stringr::str_split_i(
      stringr::str_split_i(tile_names[j], pattern = '_', i = 5), pattern = '\\.',1)
    fft_name <- paste('nx',nx,'ny',ny,'01','tile.nc', sep = '_')
    psd <- Re(f.out_stars$spectral_power_matrix) #take real part, complex part always 0i
    psd_raster <- raster::raster(psd,
                                 xmn = min(f.out$kx), xmx = max(f.out$kx),
                                 ymn = min(f.out$ky), ymx = max(f.out$ky))
    raster::writeRaster(psd_raster,
                        paste(out.dir,fft_name, sep = '/'), overwrite = T)
    p_star <- read_stars(paste(out.dir,fft_name, sep = '/'), proxy = T)
    names(p_star) <- 'power_spectrum'
    
    lst <- list()
    lst$star <- p_star
    # kx, ky dimensions can be recreated using st_dim_to_attr()
    assign(fft_name, lst, envir = fft_env)

  } else{ # split in two DEMs and generate 2 FFTs
    no_split <- no_split + 1
    n_split <- floor(dim(tile)[2]/2)
    dem_mat <- as.matrix(tile)[,1:n_split]
    f.out <- fft2D(dem_mat, dr = res(tile)[1], Hann = T)
    f.out_stars <- st_as_stars(f.out[c(1:2,5:6)])
    # in stars, rows are mapped to the first dimension, cols to the second
    names(st_dimensions(f.out_stars)) <- c('ky','kx') # kx <-> cols, consistent with fft2d
    f.out_stars <- st_set_dimensions(f.out_stars, 1, f.out_stars$ky[,1])
    f.out_stars <- st_set_dimensions(f.out_stars, 2, f.out_stars$kx[1,])
    nx <- stringr::str_split_i(tile_names[j], pattern = '_', i = 3)
    ny <- stringr::str_split_i(
      stringr::str_split_i(tile_names[j], pattern = '_', i = 5), pattern = '\\.',1)
    fft_name <- paste('nx',nx,'ny',ny,'01','tile.tif', sep = '_')
    psd <- Re(f.out_stars$spectral_power_matrix) #take real part, complex part always 0i
    psd_raster <- raster::raster(psd,
                                   xmn = min(f.out$kx), xmx = max(f.out$kx),
                                   ymn = min(f.out$ky), ymx = max(f.out$ky))
    raster::writeRaster(psd_raster,
                paste(out.dir,fft_name, sep = '/'), overwrite = T)
    p_star <- read_stars(paste(out.dir,fft_name, sep = '/'), proxy = T)
    names(p_star) <- 'power_spectrum'
    
    lst <- list()
    lst$star <- p_star
    # kx, ky dimensions can be recreated using st_dim_to_attr()
    assign(fft_name, lst, envir = fft_env)

    # repeat for second split
    dem_mat <- as.matrix(tile)[,(n_split+1):dim(tile)[2]]
    f.out <- fft2D(dem_mat, dr = res(tile)[1], Hann = T)
    fft_name <- paste('nx',nx,'ny',ny,'02','tile.nc', sep = '_')
    psd <- Re(f.out_stars$spectral_power_matrix) #take real part, complex part always 0i
    psd_raster <- raster::raster(psd,
                                 xmn = min(f.out$kx), xmx = max(f.out$kx),
                                 ymn = min(f.out$ky), ymx = max(f.out$ky))
    raster::writeRaster(psd_raster,
                        paste(out.dir,fft_name, sep = '/'), overwrite = T)
    p_star <- read_stars(paste(out.dir,fft_name, sep = '/'), proxy = T)
    names(p_star) <- 'power_spectrum'
    
    lst <- list()
    lst$star <- p_star
    # kx, ky dimensions can be recreated using st_dim_to_attr()
    assign(fft_name, lst, envir = fft_env)  }
  
}

# downstream handling -----------------------------------------------------


current_tile <- st_as_stars(p_star) #reads this from proxy to memory
current_tile
saveRDS(fft_env, file = '03_derived_data/fft_env.Rds')

length(ls(envir = fft_env))

       