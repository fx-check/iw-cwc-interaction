#run aggregating_u_v_geostrphic_to_UTM prior?
# Splitting DEM into tiles based on raster cells:

# preamble ---------------------------------------------------------------
library(terra)
library(raster)
library(stringr)
library(sf)
fn <- list.files(r'(C:\Users\fbutschek\Desktop\hycom\2022_south\S_T_data)', pattern = '[0-9].nc$',full.names =T)[1]
# nc_grid <- raster(fn)
fn <- '02_secondary_data/Blake_south_exple-20221231_21.nc'
nc_grid <-  raster(fn)
extent(nc_grid) <- c(280,284,27,29.96) #corrects shift
# getwd()
full_dem_fn <- '02_secondary_data/Bathymetry_BlakePlateau.tif'
# dem_fn <- r'(C:\Users\fbutschek\Desktop\hycom\nc_overlap_bathy_south.tif)'
(dem <- rast(full_dem_fn))
dem[dem==0]<- -NA;NAflag(dem) <- -9999
# dem <- rast(dem)
image(dem)
#reproject to Blake and crop here!

terra_grid <- terra::rast(fn, lyrs= 1)
  # nc_grid$Water.Temperature)
(nc_grid_utm <- terra::project(terra_grid,"EPSG:32618", 
                               res = c(7450, 4380),
                               method = 'bilinear')
  )
terra::compareGeom(nc_grid_utm, dem,ext=F,rowcol=F)
south_dem <- crop(dem, nc_grid_utm)
# this is a round-about way of doing it, code hasn't fully been migrated from
# raster to terra
poly_grid <- rasterToPolygons(raster(nc_grid_utm), na.rm = F)
poly_grid$water_temp_depth.0 <- NULL #remove water.temp, not needed
# Create cell index for grid ----------------------------------------------
# add x and y indices to columns and rows respectively:
for (i in 1:nrow(poly_grid@data)){ # 
  nx <- i%%dim(nc_grid_utm)[2]
  poly_grid@data$nx[i] <- ifelse(nx == 0,dim(nc_grid_utm)[2],nx)
  ny <- floor(i/dim(nc_grid_utm)[2]) # 0 for first row, 75 for last except last value per row
  end_of_row <- i%%dim(nc_grid_utm)[2]==0 # for when floor does not return lower value
  # at end of row: i.e. last column, compensate for floor division 
  #flip ny to correspond with nc convention:
  poly_grid@data$ny[i] <- ifelse(end_of_row,
                                 ny, 
                                 ny+1) #incr. by 1 except in last clmn
}


poly_sf_UTM <- st_as_sf(poly_grid)
sf::write_sf(poly_sf_UTM, 
             r'(03_derived_data\Blake/poly_grid_UTM_south_v3.shp)')


# if you want to run an example:
# i = 2870
# (DEM <- crop(dem, poly_grid_UTM[i, ]))

# Segment DEM into tiles, one per HYCOM cell ------------------------------

# loop this across all nc grid cells, and write to file:
out.dir <- r'(03_derived_data\Blake\Bathy_tiles\)'
existing_tiles <- list.files(out.dir)
for (i in 1:length(poly_grid)){
  nc_cell <- poly_grid[i, ]
  DEM <- try(crop(dem,nc_cell, snap = 'out'), silent = T)
  if(class(DEM) %in% 'try-error'){next} 
  no_na <- global(is.na(DEM), fun = 'sum')[[1]]
  if(no_na/ncell(DEM)>0.25) next # don't safe tile if more than 25% NAs
  fname <- paste0(out.dir,'tile_nx_',stringr::str_pad(poly_grid@data$nx[i],2, pad='0'),
                  '_ny_', stringr::str_pad(poly_grid@data$ny[i],2, pad='0'),'.tif')
  writeRaster(DEM, fname, overwrite = T)

} # the nx and ny cell number is included in the file name
# this corresponds to same cells in the utm nc files
# tests in UTM 18N show that the bathy tiles are aligned with the poly_grid_UTM


# Sanity Checks -----------------------------------------------------------


# some sanity checks to make sure we got the indexing right:
# nc_grid_UTM18N <- projectRaster(nc_grid, crs = crs(dem), method = 'ngb')# no good, increases cell count
plot(nc_grid)
lines(poly_grid[1:dim(nc_grid)[2],1]) #polygons are written row-wise along the 51 cols
lines(poly_grid[(dim(nc_grid)[2]+1):(2*dim(nc_grid)[2]),1])
lines(poly_grid[(dim(nc_grid)[2]*2+1):(3*dim(nc_grid)[2]),1])

# the indexing right? here the y axis is increasing up, should be increasing down,
# but that's just base R graphics making life difficult
plot(poly_grid@data$ny~poly_grid@data$nx, type = 'n')
text(poly_grid@data$nx[1:51],poly_grid@data$ny[1:51],
     labels = poly_grid@data$nx[1:51], cex = .5)
text(poly_grid@data$nx[seq(1,by = dim(nc_grid)[2], length.out = dim(nc_grid)[1])]
     ,poly_grid@data$ny[seq(1,by = dim(nc_grid)[2], length.out = dim(nc_grid)[1])],
     labels = poly_grid@data$ny[seq(1,by = dim(nc_grid)[2], length.out = dim(nc_grid)[1])], cex = .5)
bool <- poly_grid@data$ny==poly_grid@data$nx
text(x = poly_grid@data$nx[bool],
     y = poly_grid@data$ny[bool],
     label = poly_grid@data$ny[bool], cex = 0.5)

