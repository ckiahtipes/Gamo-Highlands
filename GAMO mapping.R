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

ETHP_tree = landcover("trees", path = "LAND/")
ETHP_grss = landcover("grassland", path = "LAND/")
ETHP_crop = landcover("cropland", path = "LAND/")

ETHP_tree = crop(ETHP_tree, ETHP_dem)

ESAtrsh_col=c(NA,"#6c9575")
ESAgrcp_col=c(NA,"#ffd731")
ESAcrop_col=c(NA,"darkred")

par(mfrow = c(1,2))

plot(0, xlim = c(32, 48), ylim = c(3, 15), pch = NA, axes = FALSE, ann = FALSE)
plot(ETHP_tree, add = TRUE, alpha = 0.8, col = ESAtrsh_col, legend = FALSE)
plot(ETHP_grss, add = TRUE, alpha = 0.8, col = ESAgrcp_col, legend = FALSE)
plot(ETHP_crop, add = TRUE, alpha = 0.8, col = ESAcrop_col, legend = FALSE)
plot(ETHP_dem, col = gray.colors(4661), main = "Ethiopia")

polygon(c(lon,rev(lon)),c(rep(lat, each = 2)), col = NA, lwd = 2, border = "forestgreen")

plot(GAMO_dem, xlim = lon, ylim = lat, col = gray.colors(3550), main = "Gamo Highlands Coring Locations")

points(clong, clats, pch = 21, bg = "goldenrod")

par(mfrow = c(1,1))
