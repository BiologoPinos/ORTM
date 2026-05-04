# GH_map_process
library(gtools)
library(ggplot2)
library(terra)
library(dplyr) 
library(ggrepel)
library(ggsn)
library(BMS)
library(sp)
library(sf)
library(raster)

# library(ggrepel)
# library(ggmap)
# library(reshape2)
# require(rgdal) 
# require(maptools)
# require(plyr)
# require(broom)
# require(rgeos)
# require(raster)
# require(leaflet)

ORlnd_sp = vect("./GIS/OR_polygon.shp", layer="OR_polygon")
# Define projections
utmz10prj = "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83"  # or EPSG 9807
llprj <-  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
ORlnd_sp <- project(ORlnd_sp, llprj)
# Convert SpatVector to SpatialPolygonsDataFrame using as() function from the sp package

ORlnd = st_as_sf(ORlnd_sp)
ggplot(data = ORlnd) + 
  geom_sf()
# Read in coastal segments
ORblk_sp = vect("./GIS/Sea_otter_hab_all.shp", layer="Sea_otter_hab_all")
ORblk_sp <- project(ORblk_sp, llprj)
ORblk_sf = st_as_sf(ORblk_sp)

ORblkpts = vect("./GIS/Sea_otter_hab_sct_pts.shp", layer="Sea_otter_hab_sct_pts")
ORblkpts = st_as_sf(ORblkpts)

# Read in grid
ORgrd_sp = vect("./GIS/Sea_otter_grid.shp", layer="Sea_otter_grid")
ORgrd_sp <- project(ORgrd_sp, llprj)
ORgrd_sf <- st_as_sf(ORgrd_sp)
ORgrd = ORgrd_sf
ORgrd$id <- ORgrd$gridID
ORgrd$Kdens_mn = pmin(15,ORgrd_sf$Kdens_mn)

# Merge blocks with values dataframe (use same approach for combining with density data)
val = as.character(unique(unique(ORblk_sf$ID)))
valF = as.factor(val)
values = data.frame(BlockID = valF,ID=val,
                    value = val)
ORblk = merge(ORblk_sf, values, by.x='ID', by.y = 'BlockID')
# rm(HGblk_sf)

c24 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
# cols = c24[sample(length(c24),nrow(values),replace = T)]
cols = rep(c24,3)
# Plot of blocks with labeled numbers (save this for plotting)
mapplot = ggplot() + 
  geom_sf(data=ORblk, aes(color=as.factor(ID)),
               alpha = 1,size = 1) +
  scale_fill_manual(values = cols, guide = FALSE) + 
  scale_color_manual(values = cols, guide = FALSE) + 
  geom_sf(data = ORlnd, aes(),
               color="wheat4", fill="cornsilk1",size = 0.1) +      
  # geom_text_repel(data=HGblkpts, aes(x=Xcoord, y=Ycoord, label=BlockID), 
  #                 force=4, point.padding = 0.1, box.padding = 0.5,
  #                 min.segment.length = .05,
  #                 nudge_x = -2, nudge_y = 1) +
  #scale_x_continuous(name = "Longitude", breaks = seq(-125,-123)) + # , breaks = NULL, labels = NULL
  #scale_y_continuous(name = "Latitude") + # , breaks = NULL, labels = NULL
  #north(HGlnd,location = "topright") +
  # scalebar(HGlnd, dist = 50, dist_unit = "km", st.size = 3.5, 
  #          transform = FALSE, location = "bottomleft") +  
  ggtitle("Haida Gwaii, coastal sections for sea otter model") +
  coord_sf() +
  # coord_fixed(ratio=1.5) + 
  theme_minimal()
print(mapplot)
# Example of plotting continuous values (e.g. replace value with densities)
ggplot() + 
  geom_sf(data=ORblk_sf, aes(fill = K_value, color=K_value),
               alpha = 1,size = 1) +
  scale_fill_continuous(guide = FALSE, low = "#fff7ec", high = "#7F0000") + 
  scale_color_continuous(guide = FALSE, low = "#fff7ec", high = "#7F0000") + 
  scale_x_continuous(name = "Longitude", breaks = seq(-125,-123)) + # , breaks = NULL, labels = NULL
  scale_y_continuous(name = "Latitude") + # , breaks = NULL, labels = NULL
  geom_sf(data = ORlnd, aes(fill=piece),
               color="wheat4", fill="cornsilk1",size = 0.1) +   
  # north(HGlnd,location = "topright") +
  # scalebar(HGlnd, dist = 50, dist_unit = "km", st.size = 3.5, 
  #          transform = FALSE, location = "bottomleft") +
  coord_sf() +
  # coord_equal(ratio=1) + 
  theme_minimal()
ggplot() + 
  geom_sf(data=ORgrd_sf, aes(fill = Kdens_mn, color=Kdens_mn),
               alpha = 1,size = 1) +
  scale_fill_continuous(guide = FALSE, low = "#fff7ec", high = "#7F0000") + 
  scale_color_continuous(guide = FALSE, low = "#fff7ec", high = "#7F0000") + 
  scale_x_continuous(name = "Longitude", breaks = seq(-125,-123)) + # , breaks = NULL, labels = NULL
  scale_y_continuous(name = "Latitude") + # , breaks = NULL, labels = NULL
  geom_sf(data = ORlnd, aes(fill=piece),
               color="wheat4", fill="cornsilk1",size = 0.1) +   
  # north(HGlnd,location = "topright") +
  # scalebar(HGlnd, dist = 50, dist_unit = "km", st.size = 3.5, 
  #          transform = FALSE, location = "bottomleft") +
  coord_sf() +
  # coord_equal(ratio=1) + 
  theme_minimal()
# Save mapplot and HGlnd and HGblk and HGblkpts
setwd("./data")
save(mapplot,ORlnd,ORblk,ORgrd,ORblkpts,file="GISdata.rdata")
setwd("..")


