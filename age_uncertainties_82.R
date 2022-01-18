## timing of 8.2 ka event with 8.2 ka event excursions

setwd(".../Holocene_abrupt_events/")

library(RMySQL)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

#load 8.2 ka data
signals_82 <- read.csv("spel_82_signals.csv")

# connect to SISAL database
mydb <- dbConnect(MySQL(), user = "root", password = "", dbname = "sisal_v2", 
                  host = "localhost")

# load Holocene data
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING (site_id) JOIN sample USING (entity_id) JOIN original_chronology USING (sample_id) JOIN sisal_chronology USING (sample_id) JOIN d18O USING (sample_id);")
Raw_Data <- Raw_Data %>% filter(entity_id %in% signals_82$entity_id)

dat_res <- data.frame()
noSISAL <- c()
for (i in na.omit(unique(signals_82$entity_id))){
  subdat <- signals_82 %>% filter(entity_id == i)
  
  if (nrow(Raw_Data %>% filter(sample_id == subdat$sample_id_min_bp))==0) {
    noSISAL <- c(noSISAL, i)
  } else {
    sub_start <- Raw_Data %>% filter(sample_id == subdat$sample_id_max_bp); sub_start$grp <- "start"
    sub_end <- Raw_Data %>% filter(sample_id == subdat$sample_id_min_bp); sub_end$grp <- "end"
    
    subdat2 <- rbind(sub_start, sub_end)
    subdat2 <- subdat2[,-c(1,7:10,12:36,40:42,64:66)]
    
    subdat_res <- na.omit(subdat2 %>% gather(key = "age_model", value = "age", 7:30))
    subdat_res$grp2 <- "age"
    subdat_res$grp2[which(grepl("uncert_pos", subdat_res$age_model))] <- "uncert_pos"
    subdat_res$grp2[which(grepl("uncert_neg", subdat_res$age_model))] <- "uncert_neg"
    
    subdat_res$age_model <- gsub("_age", "", subdat_res$age_model)
    subdat_res$age_model <- gsub("_uncert_pos", "", subdat_res$age_model)
    subdat_res$age_model <- gsub("_uncert_neg", "", subdat_res$age_model)
    
    subdat_res <- subdat_res %>% spread(key = grp2, value = age)
    subdat_res$uncert_neg <- subdat_res$age - subdat_res$uncert_neg
    subdat_res$uncert_pos <- subdat_res$age + subdat_res$uncert_pos
    subdat_res <- subdat_res %>% gather(key = "grp2", value = "age", 9:11)
    
    subdat_res <- subdat_res %>% spread(key = grp, value = age)
    
    dat_res <- rbind(dat_res, subdat_res)
  }
}

dat_res[which(dat_res$age_model == "interp"), "age_model"] <- "original"

##orig chron only
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING (site_id) JOIN sample USING (entity_id) JOIN d18O USING (sample_id);")
Raw_Data <- Raw_Data %>% filter(entity_id %in% noSISAL)

noSISAL2 <- signals_82 %>% filter(entity_id %in% noSISAL)

noSISAL2 <- data.frame(noSISAL2[,c("entity_id","site_id","site_name","latitude","longitude","entity_name")], 
                       age_model = "original", grp = "age", noSISAL2[,c("min_bp","max_bp")])
colnames(noSISAL2)[9:10] <- c("end","start")
colnames(dat_res)[8] <- "grp"

y <- rbind(dat_res, noSISAL2)

## regions
y$region <- with(y, ifelse(latitude >= 20 & longitude >= -10 & longitude <= 50, "Europe/Mediterranean", ifelse(
  latitude >= -10 & longitude >= 50, "Asia", ifelse(
    latitude <= 0 & longitude <= -50 & longitude <= -20, "South America", ifelse(
      latitude >0 & longitude <= -50 & longitude <= -20, "North America", "South America"
    )
  )
)))

y$site_name2 <- with(y, paste(site_name, " (", entity_name, ")", sep = ""))


## plot
p1 <- ggplot(data = filter(y, region == "Asia"), mapping = aes(xmin = end, xmax = start, y = site_name2, col = age_model, lty = grp)) + 
  geom_errorbarh(position = position_dodge()) +
  geom_vline(xintercept = 8150, lty = 2) +
  geom_vline(xintercept = c(8090,8245), lty = 3, ) +
  theme_bw() +
  ylab("") +
  expand_limits(x = c(6400,12000)) +
  theme(axis.title = element_blank(),
        legend.position = "none")
p2 <- ggplot(data = filter(y, region == "Europe/Mediterranean"), mapping = aes(xmin = end, xmax = start, y = site_name2, col = age_model, lty = grp)) + 
  geom_errorbarh(position = position_dodge()) +
  geom_vline(xintercept = 8150, lty = 2) +
  geom_vline(xintercept = c(8090,8245), lty = 3, ) +
  theme_bw() +
  ylab("")  +
  expand_limits(x = c(6400,12000)) +
  theme(axis.title = element_blank(),
        legend.position = "none")
p3 <- ggplot(data = filter(y, region == "North America"), mapping = aes(xmin = end, xmax = start, y = site_name2, col = age_model, lty = grp)) + 
  geom_errorbarh(position = position_dodge()) +
  geom_vline(xintercept = 8150, lty = 2) +
  geom_vline(xintercept = c(8090,8245), lty = 3, ) +
  theme_bw() +
  ylab("")  +
  expand_limits(x = c(6400,12000)) +
  theme(axis.title = element_blank(),
        legend.position = "none")
p4 <- ggplot(data = filter(y, region == "South America"), mapping = aes(xmin = end, xmax = start, y = site_name2, col = age_model, lty = grp)) + 
  geom_errorbarh(position = position_dodge()) +
  geom_vline(xintercept = 8150, lty = 2) +
  geom_vline(xintercept = c(8090,8245), lty = 3, ) +
  theme_bw() +
  ylab("")  +
  expand_limits(x = c(6400,12000)) +
  theme(axis.title = element_blank(),
        legend.position = "none")

pdf("age_uncert_fig3.pdf", width = 3.93, height = 9.83)
plot_grid(p2,p1,p4,p3, ncol = 1, align = "v", rel_heights = c(3.57,3.75,1.06,1.29))
dev.off()