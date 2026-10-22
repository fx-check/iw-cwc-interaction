# Load DEM tiles from grid and perform FFT with Hann filter
# Write the FFT output including wavenumber vector matrix (for subsequent masking)
# to a dedicated environment with nx ny index. This allows relating the data with
# .nc u&v values for rotation, as well as bottom N2 for energy calculation

source('04_data_analysis/02_functions/fft2d_mod.R')
library(raster)
library(stringr)
library(ncdf4)
fft_env <- new.env()
tile_dir <- "C:\\Users\\fbutschek\\Desktop\\hycom\\Bathy_tiles"
# tile_dir <- "03_derived_data/Blake/Bathy_tiles/"
tile_names <- list.files(tile_dir, pattern = '.tif$', full.names = F)

# note on dimensions: raster rows and columns correspond to matrix dimensions
# rows are y direction, columns are x directions, this is consistent with our fft2d function
# to avoid any confusion downstream, we'll be explicit about the kx ky dimension 
# and rad/s units in writing nc files using the {ncdf4} package 


# nc output file wrapper --------------------------------------------------
fft_to_nc <- function(fft_out, nc_filename, nc_output_directory){
  #define spectral power matrix nc variable with dimensions and attributes:
  kx_dim <- ncdim_def('kx', units = 'rad/m', vals = fft_out$kx[1,],
                      create_dimvar=T, longname    = "eastward wavenumber")
  ky_dim <- ncdim_def('ky', units = 'rad/m', vals = fft_out$ky[,1],
                      create_dimvar=T, longname    = "northward wavenumber")
  pds_nc <- ncvar_def(name = 'spectral_power_matrix',
                      units = 'm^4',
                      dim = list(ky_dim, kx_dim),
                      missval = NA_real_, prec = 'float')
  nc_out <- nc_create(filename = paste(nc_output_directory,
                                       nc_filename, sep = '/'), 
                      vars = pds_nc, force_v4 = TRUE)
  on.exit(nc_close(nc_out), add = T)
  ncatt_put(nc_out, 0, "Conventions", "CF-1.8")
  ncatt_put(nc_out, "kx", "axis", "X")
  ncatt_put(nc_out, "ky", "axis", "Y")
  # ncatt_put(nc_out, "kx", "long_name", "eastward wavenumber")
  # ncatt_put(nc_out, "ky", "long_name", "northward wavenumber")
  ncatt_put(nc_out, "spectral_power_matrix", "long_name", "2-D power density spectrum")
  
  ncatt_put(nc_out, 0, "fft_convention", "fft shifted kx/ky; DC off +1 index from center; Hann window Wss normalization")
  
  ncvar_put(nc_out, 'spectral_power_matrix', Re(fft_out$spectral_power_matrix))
  invisible(nc_out$filename)
}

# out.dir <- '03_derived_data/Blake/fft_tiles'
out.dir <- r'(C:\Users\fbutschek\Desktop\hycom\fft_output)'

# New tile splitting and fft computation rationale ------------------------
# here, all ncols will be split in two, unless <130 (3 cases)
unique_psd_dims <- c()
for (j in 1:length(tile_names)) {
  tile <- raster(paste(tile_dir,tile_names[j],sep = '/'))
  #skip edge cases with few rows, n = 36 or few col nums (n=8)
  nr = dim(tile)[1]; nc = dim(tile)[2]
  if(dim(tile)[1]<100 | dim(tile)[2]<200) next
  nc_split <- nc%/%2
  if(nr>128){
    pad_length <- (nr%%128)%/%2
    nr_sel <- c(rep(FALSE,pad_length), rep(TRUE, 128),
                rep(FALSE,(nr-pad_length-128)))
  }else {nr_sel <- rep(TRUE,nr)}
  NAvalue(tile)<- 0
  
  dem_mat <- as.matrix(tile)[nr_sel,] # avoids raster ambiguity in fft
  # correct N-up orientation and split:
  sub_tile <- dem_mat[nrow(dem_mat):1,1:nc_split]
  #only fft with <20% NA in sub_tile split 
  if(sum(is.na(sub_tile))/length(sub_tile) < 0.20){
    #detrend and deal with remaining NAs
    det_tile <- detrend_mat(sub_tile,
                            coords = coordinates(tile)[1:(nr*nc_split),])
    f_out <- fft2D(det_tile, dr = res(tile)[1], Hann = T) #dr = dx = dy
    #name fft based on dem tile, adding '_01' for sub-tile 1:
    nx <- str_split_i(tile_names[j], pattern = '_', i = 3)
    ny <- str_split_i(
      str_split_i(tile_names[j], pattern = '_', i = 5), 
      pattern = '\\.',1)
    tile_number <- '01'
    tile_nm <- paste('nx',nx,'ny',ny, tile_number, sep = '_')
    fft_name <- paste(tile_nm,'tile.nc', sep = '_') #replaced tiff with NC
    #write out:
    fft_to_nc(fft_out = f_out, nc_filename = fft_name,
              nc_output_directory = out.dir)
    if(paste0("tile_",tile_nm) != paste(str_split_i(tile_names[j], '.tif',1), 
                                        tile_number, sep = "_")){
      warning('file naming mismatch')}
    psd_dims <- dim(f_out$spectral_power_matrix)
    
    psd_dims <- dim(f_out$spectral_power_matrix)
    if(!(unique(psd_dims) %in% unique_psd_dims)) {
      unique_psd_dims <- c(unique_psd_dims, unique(psd_dims))}
    lst <- list('filename'= fft_name,
                'nx' = as.integer(nx),
                'ny'= as.integer(ny),
                'tile_num' = as.integer(tile_number),
                'psd_dim' = psd_dims)
    # kx, ky dimensions can be recreated

    assign(tile_nm, lst, envir = fft_env)
  }
 
  # repeat for second split
  # correct N-up orientation and split:
  sub_tile <- dem_mat[nrow(dem_mat):1,(nc_split+1):nc]
  #only fft with <20% NA in sub_tile split 
  if(sum(is.na(sub_tile))/length(sub_tile) < 0.20){
    #detrend and deal with remaining NAs
    det_tile <- detrend_mat(sub_tile,
                            coords = coordinates(tile)[(nr*nc_split+1):(nr*nc),])
    f_out <- fft2D(det_tile, dr = res(tile)[1], Hann = T) #dr = dx = dy
    #name fft based on dem tile, adding '_01' for sub-tile 1:
    nx <- str_split_i(tile_names[j], pattern = '_', i = 3)
    ny <- str_split_i(
      str_split_i(tile_names[j], pattern = '_', i = 5), 
      pattern = '\\.',1)
    tile_number <- '02'
    
    tile_nm <- paste('nx',nx,'ny',ny, tile_number, sep = '_')
    fft_name <- paste(tile_nm,'tile.nc', sep = '_')

    fft_to_nc(fft_out = f_out, nc_filename = fft_name,
              nc_output_directory = out.dir)
    
    if(paste0("tile_",tile_nm) != paste(str_split_i(tile_names[j], '.tif',1), 
                                        tile_number, sep = "_")){
      warning('file naming mismatch')}
    psd_dims <- dim(f_out$spectral_power_matrix)
    
  if(!(unique(psd_dims) %in% unique_psd_dims)) {
    unique_psd_dims <- c(unique_psd_dims, unique(psd_dims))}
  lst <- list('filename'= fft_name,
              'nx' = as.integer(nx),
              'ny'= as.integer(ny),
              'tile_num' = as.integer(tile_number),
              'psd_dim' = psd_dims)
  # kx, ky dimensions can be recreated using st_dim_to_attr()
  assign(tile_nm, lst, envir = fft_env)
  }
}


# Test run ----------------------------------------------------------------
# tiles have more columns (nc) than rows (nr) due to unequal nc spacing and
# longitude effects
(tile <- raster(paste(tile_dir,tile_names[1200],sep = '/')))
ndif <- diff(dim(tile))[1] #difference between tile dimensions
ndiff_split <- dim(tile)[1] - dim(tile)[2]/2 
#splitting makes the second dimension shorter (nc<nr), let's check which discrepancy is greater
if( ndif/min(dim(tile)[1:2]) <= ndiff_split/ max(dim(tile)[1:2])/2){ 
  #just calculate one FFT if unsplit diff smaller or equal
  dem_mat <- as.matrix(tile)
  f.out <- fft2D(dem_mat, dr = res(tile)[1], Hann = T)
  nx <- str_split_i(tile_names[1], pattern = '_', i = 3)
  ny <- str_split_i(str_split_i(tile_names[1], 
                                pattern = '_', i = 5), pattern = '\\.',1)
  fft_name <- paste('nx',nx,'ny',ny,'01', sep = '_')
  fft_to_nc(fft_out = f_out, nc_filename = fft_name, 
            nc_output_directory = out.dir)
} else{ # split in two DEMs and generate 2 FFTs
  n_split <- floor(dim(tile)[2]/2)
  dem_mat <- as.matrix(tile)[,1:n_split]
  f_out <- fft2D(dem_mat, dr = res(tile)[1], Hann = T)
  nx <- stringr::str_split_i(tile_names[1], pattern = '_', i = 3)
  ny <- stringr::str_split_i(
    stringr::str_split_i(tile_names[1], pattern = '_', i = 5), pattern = '\\.',1)
  fft_name <- paste('nx',nx,'ny',ny,'01', sep = '_')
  #spectral power nc dim and var definition:
  
  
  assign(fft_name, f.out, envir = fft_env)
  #repeat for second split
  dem_mat <- as.matrix(tile)[,(n_split+1):dim(tile)[2]]
  f_out <- fft2D(dem_mat, dr = res(tile)[1], Hann = T)
  fft_name <- paste('nx',nx,'ny',ny,'02', sep = '_')
  fft_to_nc(fft_out = f_out, nc_filename = fft_name, 
            nc_output_directory = out.dir)
  assign(fft_name, f.out, envir = fft_env)
}

# Loop across all tiles ---------------------------------------------------


#### 1.) Improved heuristics on splits:
# Cells added by padding to a square power-of-two canvas for a (nr x nc) matrix
padded_cells <- function(nr, nc) {
  L <- 2^ceiling(log2(max(nr, nc)))
  L*L - nr*nc
}

# Decide to split columns only if it reduces padding by at least `min_gain` fraction
should_split_cols <- function(nr, nc, min_gain = 0.05) {
  pad_unsplit <- padded_cells(nr, nc)
  
  nc1 <- nc %/% 2
  nc2 <- nc - nc1
  nc3 <- max(nc1,nc2)
  pad_split <- padded_cells(nr, nc3) 
  
  # Split only if it reduces padding by more than the minimum threshold set in min_gain
  gain <- (pad_unsplit - pad_split) / max(pad_unsplit, 1L)
  gain > min_gain
}
no_split = no_unsplit = 0 # record keeping of unsplit vs split matrices
unique_psd_dims <- c() #this will make downstream kx ky handling easier
# allowing to precompute identical ky,kx and facilitate flow alignment

for (j in 1:10
     # length(tile_names)
     ) {
  tile <- raster(paste(tile_dir,tile_names[j],sep = '/'))
  ndif <- diff(dim(tile))[1] #difference between tile dimensions
  ndiff_split <- dim(tile)[1] - dim(tile)[2]/2 
  
  if(!should_split_cols(dim(tile)[1], dim(tile)[2], min_gain = 0.05)
    # ndif/min(dim(tile)[1:2]) <= ndiff_split/ max(dim(tile)[1:2])/2
    ){ 
    #just calculate one FFT if unsplit diff smaller or equal
    no_unsplit <- no_unsplit+1
    dem_mat <- as.matrix(tile)
    f_out <- fft2D(dem_mat, dr = res(tile)[1], Hann = T)
    nx <- stringr::str_split_i(tile_names[j], pattern = '_', i = 3)
    ny <- stringr::str_split_i(
      stringr::str_split_i(tile_names[j], pattern = '_', i = 5), pattern = '\\.',1)
    fft_name <- paste('nx',nx,'ny',ny,'01','tile.nc', sep = '_')
    fft_to_nc(fft_out = f_out, nc_filename = fft_name, 
              nc_output_directory = out.dir)
    tile_number <- str_split_i(fft_name, '_', 5)
    tile_nm <- paste(str_split_i(tile_names[j], '.tif',1), 
                     tile_number, sep = "_")
    psd_dims <- dim(f_out$spectral_power_matrix)
    if(!(unique(psd_dims) %in% unique_psd_dims)) {
      unique_psd_dims <- c(unique_psd_dims, unique(psd_dims))}
    lst <- list('filename'= fft_name,
                'nx' = nx,
                'ny'= ny,
                'tile_num' = tile_number[i],
                'psd_dim' = psd_dims)
    # kx, ky dimensions can be recreated using st_dim_to_attr()
    assign(tile_nm, lst, envir = fft_env)

  } else{ # split in two DEMs and generate 2 FFTs
    no_split <- no_split + 1
    n_split <- floor(dim(tile)[2]/2)
    dem_mat <- as.matrix(tile)[,1:n_split]
    f_out <- fft2D(dem_mat, dr = res(tile)[1], Hann = T)
 
    nx <- stringr::str_split_i(tile_names[j], pattern = '_', i = 3)
    ny <- stringr::str_split_i(
      stringr::str_split_i(tile_names[j], pattern = '_', i = 5), pattern = '\\.',1)
    fft_name <- paste('nx',nx,'ny',ny,'01','tile.nc', sep = '_') #replaced tiff with NC
    fft_to_nc(fft_out = f_out, nc_filename = fft_name, 
              nc_output_directory = out.dir)
    tile_number <- str_split_i(fft_name, '_', 5)
    
    psd_dims <- dim(f_out$spectral_power_matrix)
    if(!(unique(psd_dims) %in% unique_psd_dims)) {
      unique_psd_dims <- c(unique_psd_dims, unique(psd_dims))}
    lst <- list('filename'= fft_name,
                'nx' = nx,
                'ny'= ny,
                'tile_num' = tile_number[i],
                'psd_dim' = psd_dims)
    # kx, ky dimensions can be recreated using st_dim_to_attr()
    assign(tile_nm, lst, envir = fft_env)

    # repeat for second split
    dem_mat <- as.matrix(tile)[,(n_split+1):dim(tile)[2]]
    f_out <- fft2D(dem_mat, dr = res(tile)[1], Hann = T)
    fft_name <- paste('nx',nx,'ny',ny,'02','tile.nc', sep = '_')
    fft_to_nc(fft_out = f_out, nc_filename = fft_name, 
              nc_output_directory = out.dir)
    tile_number <- str_split_i(fft_name, '_', 5)
    tile_nm <- paste(str_split_i(tile_names[j], '.tif',1), 
                     tile_number, sep = "_")
    psd_dims <- dim(f_out$spectral_power_matrix)
    if(!(unique(psd_dims) %in% unique_psd_dims)) {
      unique_psd_dims <- c(unique_psd_dims, unique(psd_dims))}
    lst <- list('filename'= fft_name,
                'nx' = nx,
                'ny'= ny,
                'tile_num' = tile_number[i],
                'psd_dim' = psd_dims)
    # kx, ky dimensions can be recreated using st_dim_to_attr()
    assign(tile_nm, lst, envir = fft_env)
  }

}



# downstream handling -----------------------------------------------------

# saveRDS(fft_env, file = '03_derived_data/fft_env.Rds')
# length(ls(envir = fft_env))


# recover fft_env from filenames --------------------------------------
fft_env<- new.env()

OD <- "C:/Users/fbutschek/OneDrive - University College Cork/"
OD_subdir <- "CullenPhD/16_NEA_coral_provinces/19_git/iw-cwc-interaction"
setwd(paste0(OD, OD_subdir))

all_ffts <- list.files('03_derived_data/Blake/fft_tiles/')
unique_psd_dims <- 128
for (i in all_ffts){
  a <- str_split_i(i,'\\.',1)
  b <- str_split_1(a,'_')
  nx <- as.integer(b[2]); ny = as.integer(b[4]); tile_number = as.integer(b[5])
  psd_dims <- ifelse(length(unique_psd_dims)==1, unique_psd_dims, NA)
  tile_nm <-   str_split_i(a,'_tile',1)
  lst <- list('filename'= i,
              'nx' = nx,
              'ny'= ny,
              'tile_num' = tile_number,
              'psd_dim' = psd_dims)
  assign(x = tile_nm, value = lst, envir = fft_env)
}

# Add grid convergence offset to fft_env ----------------------------------


gamma_nc <- nc_open(r'(C:\Users\fbutschek\Desktop\hycom\2022_south\gamma_grid.nc)')
gamma_grid <- ncvar_get(gamma_nc, 'gamma')
# all(dplyr::near(t(gamma_grid2)-gamma_grid, 0, tol = 1e-6)) #earlier testing
nc_close(gamma_nc)
for (i in ls(fft_env)){
  lst <- get(i, fft_env)
  lst$gamma<- gamma_grid[lst$nx, lst$ny]
  assign(i, lst, fft_env)
}
