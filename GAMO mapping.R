###Gamo Mapping

#Libraries

library(geodata)
library(maps)
library(grDevices)
library(sf)
library(raster)
library(rnaturalearth)

#Cores

cores = c("CHA", "CHO", "KAO")
clats = c(6.465545, 6.388072, 6.416460)
clong = c(37.643475, 37.668716, 37.633934)


#Simple DEM test

lat = c(6.3, 6.6)
lon = c(37.5, 37.8)

if(!file.exists("SRTM/elevation/ETH_elv_msk.tif")){
  ETHP_dem <- elevation_30s(country="ETHIOPIA", path = "SRTM/")
} else {
  ETHP_pth <- "SRTM/elevation/ETH_elv_msk.tif"
  ETHP_dem <- rast(ETHP_pth) 
}

if(!file.exists("SRTM/elevation/srtm_44_11.tif")){
  GAMO_dem <- elevation_3s(lon = 35, lat = 10, path = "SRTM/")
} else {
  GAMO_pth <- "SRTM/elevation/srtm_44_11.tif"
  GAMO_dem <- rast(GAMO_pth) 
}

#Let's use the country map to extract data more precisely.

ETHP_box = st_bbox(ETHP_dem)

#Let's pull some of the landcover data

#Trees!

if(!file.exists("LAND/landuse/WorldCover_trees_30s.tif")){
  ETHP_tree = landcover("trees", path = "LAND/")
} else {
  tree_path <- "LAND/landuse/WorldCover_trees_30s.tif"
  ETHP_tree <- rast(tree_path) 
}

#Grassland!

if(!file.exists("LAND/landuse/WorldCover_grassland_30s.tif")){
  ETHP_grss = landcover("grassland", path = "LAND/")
} else {
  grass_path <- "LAND/landuse/WorldCover_grassland_30s.tif"
  ETHP_grss <- rast(grass_path) 
}

#Cropland!

if(!file.exists("LAND/landuse/WorldCover_cropland_30s.tif")){
  ETHP_crop = landcover("cropland", path = "LAND/")
} else {
  crop_path <- "LAND/landuse/WorldCover_cropland_30s.tif"
  ETHP_crop <- rast(crop_path) 
}

#Crop around Ethopia

ETHP_tree = crop(ETHP_tree, ETHP_dem, mask = TRUE)
ETHP_tree[ETHP_tree < 0.1] <- NA

ETHP_grss = crop(ETHP_grss, ETHP_dem, mask = TRUE)
ETHP_grss[ETHP_grss < 0.1] <- NA

ETHP_crop = crop(ETHP_crop, ETHP_dem, mask = TRUE)
ETHP_crop[ETHP_crop < 0.1] <- NA

ESAtrsh_col=c("#6c9575")
ESAgrcp_col=c("#ffd731")
ESAcrop_col=c("darkred")

par(mfrow = c(1,2))

plot(ETHP_dem, col = gray.colors(4661), main = "Ethiopia", legend = FALSE)
plot(ETHP_tree, add = TRUE, col = ESAtrsh_col, legend = FALSE, alpha = 0.5)
plot(ETHP_grss, add = TRUE, col = ESAgrcp_col, legend = FALSE, alpha = 0.5)
plot(ETHP_crop, add = TRUE, col = ESAcrop_col, legend = FALSE, alpha = 0.5)


polygon(c(lon,rev(lon)),c(rep(lat, each = 2)), col = NA, lwd = 2, border = "black")

plot(GAMO_dem, xlim = lon, ylim = lat, col = gray.colors(3550), main = "Gamo Highlands Coring Locations")
plot(ETHP_tree, add = TRUE, col = ESAtrsh_col, legend = FALSE, alpha = 0.5)
plot(ETHP_grss, add = TRUE, col = ESAgrcp_col, legend = FALSE, alpha = 0.5)
plot(ETHP_crop, add = TRUE, col = ESAcrop_col, legend = FALSE, alpha = 0.5)

points(clong, clats, pch = 21, bg = "goldenrod")

par(mfrow = c(1,1))
