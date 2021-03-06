## plot fig 1: map

setwd(".../spel_abrupt_Holocene")

library(rgdal)
library(ggplot2)
library(dplyr)
library(ggnewscale)


# load world map data (downloaded from https://www.naturalearthdata.com/downloads/110m-physical-vectors/)
wmap <- readOGR(dsn = "ne_110m_land", layer = "ne_110m_land")
wmap@data$id <- rownames(wmap@data)
worldMap <- fortify(wmap)
wmap_DF <- merge(worldMap, wmap@data, by = "id")


# load WOKAM (carbonate bedrock) data (downloaded from https://produktcenter.bgr.de/terraCatalog/OpenSearch.do?search=473d851c-4694-4050-a37f-ee421170eca8&type=/Query/OpenSearch.do)
wokam <- readOGR(dsn = "WHYMAP_WOKAM/shp", layer = "whymap_karst__v1_poly")
wokam@data$id <- rownames(wokam@data)
wokamMap <- fortify(wokam)
wokam_DF <- merge(wokamMap, wokam@data, by = "id")


# load spel sites
sites_82 <- read.csv("spel_82_signals.csv")
sites_82_nosignal <- read.csv("spel_nosignal_82.csv")
sites_82 <- rbind(sites_82[,1:6], sites_82_nosignal)
sites_Hol <- read.csv("C:/Users/ph805612//OneDrive - University of Reading/Documents/abrupt_Holocene/Hol_bp_nentities.csv")
sites_Hol <- sites_Hol %>% group_by(site_id, entity_id, longitude, latitude) %>% summarise(n()) 

sites_Hol <- sites_Hol %>% filter(!site_id %in% sites_82$site_id)

sites_82$grp <- "8.2 ka analysis"
sites_Hol$grp <- "Holocene analysis"

all_sites <- rbind(sites_Hol[,c(1:4,6)], sites_82[,c(1,3,5:7)])

all_sites$grp <- factor(all_sites$grp, c("Holocene analysis", "8.2 ka analysis"))

pdf(file = "map1_fig.pdf", width = 15/2.54, height = 7/2.54)
ggplot() +
  geom_polygon(data = wmap_DF, aes(x = long, y = lat, group = group), fill = "#eeeeee") +
  
  geom_polygon(data = wokam_DF, aes(x = long, y = lat, fill = "Carbonate bedrock", group = id)) +
  scale_fill_manual(values = "#FDAE6B") + ## find paler shade of orange
  
  geom_path(data = wmap_DF, aes(x = long, y = lat, group = group), col = "dark grey", size = 0.2) +
  
  new_scale_fill() +
  geom_point(data = all_sites, aes(x = longitude, y = latitude, fill = factor(grp, levels = c("Holocene analysis", "8.2 ka analysis")), shape = factor(grp, levels = c("Holocene analysis", "8.2 ka analysis"))), size = 1.5, stroke = 0.2) +
  scale_fill_manual(values = c("#5D3A9B","#73F1C5")) +
  scale_shape_manual(values = c(21,24)) +
  
  theme_bw() +
  coord_cartesian(ylim = c(-50,65), xlim = c(-150,170)) +
  theme(axis.title = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "#F6FFFF"),
        axis.text = element_blank(),
        axis.ticks.length = unit(-0.1, "cm"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.position = c(0.13, 0.25))
dev.off()
