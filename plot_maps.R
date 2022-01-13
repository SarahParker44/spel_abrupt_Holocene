### Plot 8.2 ka signals

setwd(".../speleothem_8_2_kyr_signals")

library(rgdal)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

## load world map data (downloaded from https://www.naturalearthdata.com/downloads/110m-physical-vectors/)
wmap <- readOGR(dsn = "ne_110m_land", layer = "ne_110m_land")
wmap@data$id <- rownames(wmap@data)
worldMap <- fortify(wmap)
wmap_DF <- merge(worldMap, wmap@data, by = "id")


## plot 1: speleothem d18Oanomalies

# load spel d18O anomalies from breakpoint analysis
spel_82 <- read.csv("spel_82_signals.csv")
spel_82 <- spel_82[-which(spel_82$entity_id == 449),]
spel_no82 <- read.csv("spel_nosignal_82.csv")

# load spel anomalies from literature (data not available in SISAL or online)
non_avail_spel <- read.csv("spel_82_signals_notavailable.csv")

# plot
p1 <- ggplot() +
  
  geom_polygon(data = wmap_DF, aes(x = long, y = lat, group = id), fill = "light grey") +
  
  geom_point(data = spel_no82, aes(x = longitude, y = latitude, shape = "no change"), size = 2, col = "#919191") +
  scale_shape_manual(values = 16) +
  
  geom_point(data = spel_82, aes(x = longitude, y = latitude, fill = anom), shape = 21, size = 3) +
  scale_fill_fermenter(palette = "PiYG", n.breaks = 10, limits = c(-1,1)) +
  
  geom_point(data = non_avail_spel, aes(x = longitude, y = latitude, col = signal), shape = 4, stroke = 2) +
  scale_color_manual(values = c("#4D9221","#C51B7D")) +
  
  coord_cartesian(xlim = c(-170,180), ylim = c(-55,78)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(-180,180,30), minor_breaks = seq(-180,180,10)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(-90,90,30), minor_breaks = seq(-90,90,10)) +
  
  theme_bw() +
  theme(legend.position = c(0.1,0.3),
        legend.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank())

pdf("spel_anom_map.pdf", width = 25/2.54, height = 12/2.54)
p1
dev.off()

## plot 2: plot speleothem (non-d18O) signals

# load signals
non_d18O <- read.csv("other_spel_evidence.csv")
non_d18O_nochange <- non_d18O %>% filter(Interp == "no change")
non_d18O <- non_d18O %>% filter(Interp != "no change")
hiatuses <- read.csv("hiatus_82_dat.csv")

# plot
p2 <- ggplot() +
  
  geom_polygon(data = wmap_DF, aes(x = long, y = lat, group = id), fill = "light grey") +
  
  #geom_line(mapping = aes(x = c(110.2500,108.5), y = c(30.2700, 25.17)), col = "#BF812D", size = 2) +
  geom_point(data = non_d18O, aes(x = longitude, y = latitude, fill = Interp),size = 3, shape = 21) +
  geom_point(data = non_d18O_nochange, aes(x = longitude, y = latitude, shape = "no change"), col = "#919191", size = 2) +
  scale_shape_manual(values = 16) +
  scale_fill_manual(values = c("#BF812D","#35978F")) +
  
  geom_text(data = hiatuses, aes(x = longitude, y = latitude, col = "hiatus"), label = "H") +
  scale_color_manual(values = "black") +
  
  coord_cartesian(xlim = c(-170,180), ylim = c(-55,78)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(-180,180,30), minor_breaks = seq(-180,180,10)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(-90,90,30), minor_breaks = seq(-90,90,10)) +
  
  theme_bw() +
  theme(legend.position = c(0.1,0.2),
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
  

  


## plot 3: Morrill et al., 2013 precip anomalies

# load anomalies
library("readxl")
precip_anoms <- read_excel("Morrill_precip.xlsx")

# reshape
no_precip_anom <- precip_anoms %>% filter(anomaly == "no change")
precip_anoms <- precip_anoms %>% filter(anomaly != "no change")

# plot
p3 <- ggplot() +
  
  geom_polygon(data = wmap_DF, aes(x = long, y = lat, group = id), fill = "light grey") +
  
  geom_point(data = no_precip_anom, mapping = aes(x = longitude, y = latitude, col = "no change"), size = 2) +
  scale_color_manual(values = "#919191") +
  geom_point(precip_anoms, mapping = aes(x = longitude, y = latitude, fill = anomaly), shape = 21, size = 3) +
  scale_fill_manual(values = c("#BF812D", "#35978F")) +
  
  coord_cartesian(xlim = c(-170,180), ylim = c(-55,78)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(-180,180,30), minor_breaks = seq(-180,180,10)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(-90,90,30), minor_breaks = seq(-90,90,10)) +
  
  theme_bw() +
  theme(legend.position = c(0.1,0.2),
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())



## plot 4: Morrill et al., 2013 temp anomalies

# load anomalies
temp_anoms <- read_excel("Morrill_temp.xlsx")

# reshape
no_temp_anom <- temp_anoms %>% filter(anomaly == "no change")
temp_anoms <- temp_anoms %>% filter(anomaly != "no change")

# plot
p4 <- ggplot() +
  
  geom_polygon(data = wmap_DF, aes(x = long, y = lat, group = id), fill = "light grey") +
  
  geom_point(data = no_temp_anom, mapping = aes(x = longitude, y = latitude, col = "no change"), size = 2) +
  scale_color_manual(values = "#919191") +
  geom_point(temp_anoms, mapping = aes(x = longitude, y = latitude, fill = anomaly), shape = 21, size = 3) +
  scale_fill_manual(values = c("#4393C3","#D6604D")) +
  
  coord_cartesian(xlim = c(-170,180), ylim = c(-55,75)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(-180,180,30), minor_breaks = seq(-180,180,10)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(-90,90,30), minor_breaks = seq(-90,90,10)) +
  
  theme_bw() +
  theme(legend.position = c(0.1,0.2),
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())


## multiplot
pdf("82_maps.pdf", width = 20/2.54, height = 10/2.54)
plot_grid(p1,p3,p2,p4, ncol = 2, align = "hv", labels = c("a)","c)","b)","d)"))
dev.off()


