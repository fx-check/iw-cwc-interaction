# Splitting DEM into tiles based on raster cells:
library(raster)
library(stringr)
fn <- list.files(r'(C:\Users\fbutschek\Desktop\hycom\2020)', pattern = '.nc$',full.names =T)[1]
nc_grid <- raster(fn)
dem_fn <- r'(C:\Users\fbutschek\Desktop\hycom\nc_overlap_bathy.tif)'
(dem <- raster(dem_fn))
NAvalue(dem) <- 0
image(dem)
poly_grid <- rasterToPolygons(nc_grid, na.rm = F)
# add x and y indices to collumns and rows respectively:
for (i in 1:nrow(poly_grid@data)){ # 
  nx <- i%%dim(nc_grid)[2]
  poly_grid@data$nx[i] <- ifelse(nx == 0,dim(nc_grid)[2],nx)
  ny <- floor(i/dim(nc_grid)[2]) + 1
  endrow <- i%%dim(nc_grid)[2]==0
  poly_grid@data$ny[i] <- ifelse(endrow,ny-1,ny)
}
poly_grid$Water.Temperature <- NULL #remove water.temp, not needed

#let's run and example:
n = 2870
poly_grid_UTM <- spTransform(poly_grid, CRSobj = crs(dem))
(DEM <- crop(dem, poly_grid_UTM[n, ]))
# loop this across all nc grid cells, and write to file:
out.dir <- r'(C:\Users\fbutschek\Desktop\hycom\Bathy_tiles/)'
for (i in 1:length(poly_grid_UTM)){
  nc_cell <- poly_grid_UTM[i, ]
  DEM <- try(crop(dem,nc_cell), silent = T)
  if(class(DEM) %in% 'try-error'){next} 
  no_na <- cellStats(is.na(DEM), stat = sum)
  if(no_na/ncell(DEM)>0.15) next # don't safe tile if more than 15% NAs
  fname <- paste0(out.dir,'tile_nx',str_pad(poly_grid_UTM@data$nx[i],2, pad='0'),
                  '_ny_', str_pad(poly_grid_UTM@data$ny[i],2, pad='0'),'.tif')
  writeRaster(DEM, fname, overwrite = T)

} # the nx and ny cell number is included into the file name


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

