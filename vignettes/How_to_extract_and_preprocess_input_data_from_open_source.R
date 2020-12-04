## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load_package_data--------------------------------------------------------
library(urbandrain)
library(osmdata)

# access provided boundary polygon from example data:
data(ex1, package = "urbandrain")
boundary_polygon <- ex1$boundary_polygon


## ----osm_data_processing, warning=FALSE, message=FALSE------------------------
# rectangular area for OSM extraction:
coords <- as.data.frame(sf::st_coordinates(boundary_polygon))
min_x <- min(coords$X)
max_x <- max(coords$X)
min_y <- min(coords$Y)
max_y <- max(coords$Y)

# function to extract streets, buildings and landuse from OSM:
extractFromOSM <- function(key, boundary_polygon, min_x, min_y, max_x, max_y){
  # download street data from OSM:
  data <- osmdata::opq(bbox = c(min_x, min_y, max_x, max_y))
  data <- osmdata::add_osm_feature(opq = data, key = key)
  data <- osmdata::osmdata_sf(data)
  
  # extract polylines or polygon data:
  if(key == "highway"){
    data <- data$osm_lines[c("highway", "maxspeed")]
  }
  if(key == "landuse" | key == "building"){
    data <- data$osm_polygons[key]
  }
  
  # cut to model area
  data <- sf::st_intersection(boundary_polygon, data)
  
  return(data)
  
}

# download data from OSM:
streets <- extractFromOSM("highway", boundary_polygon, min_x, min_y, max_x, max_y)
building <- extractFromOSM("building", boundary_polygon, min_x, min_y, max_x, max_y)
landuse <- extractFromOSM("landuse", boundary_polygon, min_x, min_y, max_x, max_y)

# select streets under which the drainage network shall be placed
# for example exclude motorways, paths etc.
streets_selected <- streets[streets$highway %in% c("unclassified","residential","tertiary",
                                                   "secondary","primary",
                                                   "pedestrian","trunk",
                                                   "unclassified", "living_street"),]
# ... plot data
plot(streets["highway"])


# union landuse and building data:
a <- landuse[,2]
b <- building[,2]
colnames(b)  <- c("landuse", "geometry")
op1 <- sf::st_difference(a,sf::st_union(b)) #notice the use of st_union()
op2 <- sf::st_intersection(b, a) #notice the order of b and a
op2 <- op2[,1]
colnames(op2) <- c("landuse", "geometry")
landuse_complete <- rbind(op1, op2)

# ... plot data
plot(landuse_complete)

# additionally a data.frame specifying the imperviousness of the landuse classes is needed, for example:
landuse_classes <- data.frame(class = unique(landuse_complete$landuse), imperviousness = c(0.4,0,1))
landuse_classes


## ----srtm_dtm_processing------------------------------------------------------
# load srtm dtm data from external data provided in the package:
dtm_path <- system.file("extdata", "srtm_germany_dtm.tif", package = "urbandrain")

# load srtm tif:
dtm_germany <- raster::raster(dtm_path)

# cut raster extent to boundary_polygon:
dtm <- raster::intersect(dtm_germany, boundary_polygon)

# if needed: interpolate to get a finer raster resolution:
dtm_fine <- raster::disaggregate(dtm, fact = 10, method = "bilinear")

# ... plot data
raster::plot(dtm_fine)


