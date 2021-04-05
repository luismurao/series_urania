#' Project: Migratory patterns of Urania boisduvalii (Lepidoptera: Uraniidae)
#' as a function of biotic interactions.
#' Code to convert the shapefiles of suitability to raster. It uses the shape
#' of Cuba.
#' The script is designed to manage memory as shapefile are large
#' Paper authors: Claudia Nunez-Penichet
#'                Jorge Soberon
#                 Luis Osorio-Olvera
#' Code author: Luis Osorio-Olvera
library(rgdal)
library(raster)

#' Read shapefiles. Note that by using sf::st_read the reading could be
#' even faster.

cuba <- rgdal::readOGR("shapefiles","Cuba_vector")
occ_ura <- rgdal::readOGR("shapefiles","ura_mil_an_ag_cob_v")
occ_omp <- rgdal::readOGR("shapefiles","omp_water_an_mil_soils_veg_pro_cov")
#' Create a raster file withe the shape of Cuba. This is used as a mask to
#' generate Urania and Omphalea models
cuba_grd <- bam::shape2Grid(cuba,0.008333333)
cuba_grd_ura <- cuba_grd*0
cuba_grd_omp <-  cuba_grd*0
# IDs of suitable cells for Urania
valcell_ura<- lapply(occ_ura$OBJECTID, function(x){
  print(x)
  C_EXT <- raster::extract(cuba_grd,
                           occ_ura[occ_ura$OBJECTID==x,],
                           cellnumbers=TRUE)
  print(C_EXT)
  return(C_EXT)
})

# Fill the mask with ones for suitable pixels
for (i in 1:length(valcell_ura)) {
  print(i)
  cuba_grd_ura[valcell_ura[[i]][[1]][,1]] <- 1

}

raster::plot(cuba_grd_ura)

writeRaster(cuba_grd_ura,"modelo_urania.tif")

# IDs of suitable cells for Omphalea
valcell_omp<- lapply(occ_omp$ID, function(x){
  print(x)
  C_EXT <- raster::extract(cuba_grd,
                           occ_omp[occ_omp$ID==x,],
                           cellnumbers=TRUE)
  print(C_EXT)
  return(C_EXT)
})


# Fill the mask with ones for suitable pixels
for (i in 1:length(valcell_omp)) {
  print(i)
  cuba_grd_omp[valcell_omp[[i]][[1]][,1]] <- 1

}

writeRaster(cuba_grd_omp,"modelo_omphalea.tif",overwrite=T)

raster::plot(cuba_grd_omp)

s1 <- raster::stack(cuba_grd_ura,cuba_grd_omp)
x11()
raster::plot(s1)
