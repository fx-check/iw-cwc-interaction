

# ===================== Depth-integrated Body Force =========================

# This script computes depth-integrated body force normalized to tidal flux Q 
# from (a) depth-integrated N2 (b) seabed depth H and (c) directional slope 
# (1/h)_x along the M2 smj axis. Directional slope h'(x) was previously smoothed  
# with a depth-aware and speed-aware Gaussian filter to ensure L>>H and L>1U\omega
# here, we smooth depth lightly for the implementation of the Baines formulation
# of slope (1/h)_x = h'(x)/h^2 using chain rule derivation
# $B = F_{bd}/Q = \frac{-sin \omega t}{\omega}(\frac{1}{h})_x zN^2$
# B_int = \frac{\alpha}{\omega H^2} \int_0^H{z N^2(z) dz}
library(terra)
library(MultiscaleDTM)
library(TideHarmonics)
# ====================== INPUT DATA AND PARAMS ==============================
out_dir <- "03_derived_data/NovaScotia/output/"
# angular frequency ω for m2 tide in rad/s 
o_m2 <- harmonics[which(harmonics$name=='M2'),'speed']/180*pi/3600
NS_bathy <- rast('clipped_bathy.tif')
cellsize_m <- res(NS_bathy)[1]
k <- .5 #50m sigma, below implemented with matrix returned to 4*sigma:
(light_smooth <- focalMat(NS_bathy, d = c(k*cellsize_m,4*k*cellsize_m),type = "Gauss")  )
bathy_smooth <- focal(NS_bathy, light_smooth, 'sum')

slope_smooth <- rast(r'(03_derived_data\\smoothed_slopes_10H/slope_smooth_blended_10H.tif)')
q_flux <- rast('Qflux_projected.tif')
## Load depth-integrated zN2  and calculate body force ---------------------
season <- c("winter"=1,
            "summer"=2,
            "annual"=3)
in_dir <- r'(03_derived_data\NovaScotia\intermediary_files)'

fnames <- c("winter"=paste(in_dir,"zN2_winter_100m_fullWC.tif",sep = "/"),
            "summer"=paste(in_dir,"zN2_summer_100m_fullWC.tif",sep = "/"),
            "annual"=paste(in_dir,"zN2_100m_fullWC.tif", sep = "/")
            )
### Create N2 mask where N<ω -------------------------------------------------
nearbedN <- sqrt(rast(paste(in_dir,"N2_nearbed_100m_Combined.tif",sep = "/")))
plot(nearbedN>o_m2 | nearbedN>0)
# the upper-frequency limit (N) does not obstruct the generation of M2 IWs
# 

winter_nearbedN2 <- raster::raster(paste(in_dir,"N2_winter_100m_20percH_200m.tif",
                               sep = "/"))
(inst_rcl <- matrix(data = c(-Inf,0,0, Inf,0,1), 
                    nrow = 2,dimnames = list(c(),
                                             c('from', 'to', 'becomes'))))
(rcl_mat <- matrix(data = c(-Inf,1e-20,o_m2,0,(o_m2-1e-20), Inf,0,1,0), nrow = 3,dimnames = list(c(),c('from', 'to', 'becomes'))))
winter_mask <- raster::reclassify(winter_nearbedN2,rcl =rcl_mat)
writeRaster(winter_mask,paste(in_dir,'winter_N2mask.tif', sep = "/"))
### Calculate body force -----------------------------------------------------
for(i in season){
  season_name <- names(season)[i]
  (zN2 <- rast(fnames[[i]]))
  (B_int <- abs(slope_smooth) / (o_m2 * bathy_smooth^2) * zN2)
  (F_bd_int <- q_flux* abs(slope_smooth) / (o_m2 * bathy_smooth^2) * zN2)
  
  names(B_int) <- paste(season_name,'depth_integrated_body_force_normalised_to_Q',
                        sep = "_")
  plot(B_int, main = gsub('_'," ", names(B_int)))
  # Raster outputs were compressed for file upload only, commented out 
  # as this will cause rounding errors
  writeRaster(B_int,paste0(out_dir,
                           paste(season_name,
                                 'depth_integrated_Bd_V2.tif', 
                                 sep = "_")),
              overwrite = T
              # , wopt = list(
              #   datatype = "FLT4S",  
              #   gdal = c(
              #     "COMPRESS=ZSTD",   
              #     "PREDICTOR=3",     
              #     "ZSTD_LEVEL=15"    
              #   )
              # )
              )
 
  
  names(F_bd_int) <- paste(season_name,'depth_integrated_body_force_J_kg-1',
                           sep = "_")
  
  writeRaster(F_bd_int, 
              paste0(out_dir,
                     paste(season_name,
                           'depth_integrated_F_bd_J_kg-1.tif',
                           sep = "_")),
              overwrite = T
              # ,wopt = list(
              #   datatype = "FLT4S",  
              #   gdal = c(
              #     "COMPRESS=ZSTD",   
              #     "PREDICTOR=3",     
              #     "ZSTD_LEVEL=15"    
              #   )
              # )
              )
}

