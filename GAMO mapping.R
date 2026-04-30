###Gamo Mapping

#Libraries

library(geodata)
library(maps)
library(grDevices)
library(sf)
library(raster)
library(rnaturalearth)

#Set logicals for saving figures

save_figs = TRUE

#Cores

settl = c("Ochollo Mullato", "Garu Shongalle", "Acoma Lasha Chilashe")
slats = c(6.460783, 6.383350, 6.433733)
slong = c(37.653417,37.676950, 37.605600)

villg = c("Gibe Mullato", "Mogesa Shonggalle", "Teka Chilashe")
vlats = c(6.460950, 6.390867, 6.410983)
vlong = c( 37.645467,  37.674267, 37.625767)

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

#Attempting topo lines

topo_breaks = seq(0, 3000, 1000)

ETHP_topo = as.contour(ETHP_dem, levels = topo_breaks)

if(save_figs == TRUE){
  setEPS()
  tiff("Figure-1_maps.tiff", height = 3000, width = 2000, res = 300)
}

par(mfrow = c(2,1))

plot(ETHP_dem, col = gray.colors(4661), main = "Ethiopia", legend = FALSE)
plot(ETHP_tree, add = TRUE, col = ESAtrsh_col, legend = FALSE, alpha = 0.5)
plot(ETHP_grss, add = TRUE, col = ESAgrcp_col, legend = FALSE, alpha = 0.5)
plot(ETHP_crop, add = TRUE, col = ESAcrop_col, legend = FALSE, alpha = 0.5)
plot(ETHP_topo, add = TRUE, lwd = 0.75, alpha = 0.5)

polygon(c(lon,rev(lon)),c(rep(lat, each = 2)), col = NA, lwd = 2, border = "black")

topo_breaks = seq(0, 3000, 200)

ETHP_topo = as.contour(ETHP_dem, levels = topo_breaks)

legend(45, 14, c("Forest", "Savanna", "Cropland"), pch = c(22,22,22), pt.bg = c("#6c9575","#ffd731","darkred"), cex = 0.5)

plot(GAMO_dem, xlim = lon, ylim = lat, col = gray.colors(3550), main = "Gamo Highlands Coring Locations", legend = FALSE)
plot(ETHP_tree, add = TRUE, col = ESAtrsh_col, legend = FALSE, alpha = 0.5)
plot(ETHP_grss, add = TRUE, col = ESAgrcp_col, legend = FALSE, alpha = 0.5)
plot(ETHP_crop, add = TRUE, col = ESAcrop_col, legend = FALSE, alpha = 0.5)
plot(ETHP_topo, add = TRUE, lwd = 0.75, alpha = 0.5)

points(clong, clats, pch = 21, bg = "skyblue", cex = 1.2)
text(clong, clats+0.0175, cores)
points(vlong, vlats, pch = 22, bg = "violet", cex = 1.2)
#text(vlong-0.05, vlats, villg)
points(slong, slats, pch = 23, bg = "white", cex = 1.2)
#text(slong, slats-0.01, settl)

legend(37.725, 6.58, c("Cores", "Settlement", "Historic Settlement"), pch = c(21, 22, 23), pt.bg = c("skyblue", "violet", "white"), cex = 0.5)

if(save_figs == TRUE){
  dev.off()
}

par(mfrow = c(1,1))
