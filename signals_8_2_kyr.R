## Constrain 8.2 kyr signals 

setwd(".../speleothem_8_2_kyr")

library(dplyr)
library(strucchange)
library(ggplot2)
library(RMySQL)

# connect to SISAL database
mydb <- dbConnect(MySQL(), user = "root", password = "", dbname = "sisal_v2", 
                  host = "localhost")

# load 7.4 to 9.0 ka  data
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING (site_id) JOIN sample USING (entity_id) JOIN original_chronology USING (sample_id) JOIN d18O USING (sample_id)
                       WHERE (interp_age BETWEEN 7400 AND 9000);")
Raw_Data <- Raw_Data %>% filter(entity_status != "superseded")

## Select entities with a sufficient resolution between 7800 and 8400 years (of 30 years)
min_res <- 30

# load function for calculating mean temporal resolution (excluding hiatuses and gaps) - 'get_ent_sampling'
source("speleothem_8_2_kyr_signals/entity_sampling_mean_res.R")

res_out <- data.frame()
for (i in unique(Raw_Data$entity_id)){ # for every entity
  subdat <- Raw_Data %>% filter(entity_id == i & interp_age >= 7800 & interp_age <= 8400)
  if (nrow(subdat) <= 1){ next }
  mean_res <- get_ent_sampling(entity_id = i, age_start = 7800, age_end = 8400)$sampling_mean
  
  res_out <- rbind(res_out,
                   data.frame(unique(subdat[c("site_id","site_name","entity_id")]),
                              mean_res = mean_res))
}

res_out <- res_out %>% filter(mean_res <= 30)

## filter to records of a sufficient resolution
res_out$record_length <- NA
for (i in 1:length(res_out$entity_id)){
  ent_id <- res_out[i,"entity_id"]
  subdat <- Raw_Data %>% filter(entity_id == ent_id & interp_age >=7800 & interp_age <= 8400)
  sub_length <- max(subdat$interp_age) - min(subdat$interp_age)
  res_out[i,"record_length"] <- sub_length
}
res_out <- res_out %>% filter(record_length >= 300)

## breakpoint analysis

dtrend_dat <- data.frame()
bp_dat <- data.frame()
no_signal <- data.frame()
# breakpoint analysis
for (i in unique(res_out$entity_id)){
  subdat <- Raw_Data %>% filter(entity_id == i & interp_age >= 7800 & interp_age <= 8400)
  
  if (nrow(subdat) <= 15){ next }
  
  ## detrend
  subdat_lm <- lm(d18O_measurement ~ interp_age, data = subdat)
  lm_predicted <- predict(subdat_lm)
  subdat$detrended_d18O <- residuals(subdat_lm)
  
  bp <- breakpoints(subdat$detrended_d18O ~ 1) #breakpoint analysis

  if (length(bp$breakpoints) == 0){ 
    no_signal <- rbind(no_signal, unique(subdat[,c("site_id","site_name","entity_id","entity_name","longitude","latitude")])) }
  else if (is.na(length(bp$breakpoints))) { opt_brks <- as.numeric(names(which.min(bpts_sum$RSS["BIC",]))) } 
  else {
    ci_x <- confint(bp, breaks = length(bp$breakpoints)) #get timings of bp's with conf intervals
    
    if (any(ci_x$confint[,c(1,3)] == 0) | any(ci_x$confint[,c(1,3)] < 0) | any(is.na(ci_x$confint[,c(1,3)]))){ 
      no_signal <- rbind(no_signal, unique(subdat[,c("site_id","site_name","entity_id","entity_name","longitude","latitude")])) 
    } else {
      # output breakpoints
      ci_ages <- data.frame(site_id = unique(subdat$site_id), site_name = unique(subdat$site_name), entity_id = unique(subdat$entity_id), entity_name = unique(subdat$entity_name),
                            longitude = unique(subdat$longitude), latitude = unique(subdat$latitude),
                            bp = subdat$interp_age[ci_x$confint[,2]],
                            CI2_5 = subdat$interp_age[ci_x$confint[,1]],
                            CI97_5 = subdat$interp_age[ci_x$confint[,3]]
      )
      
      bp_dat <- rbind(bp_dat, ci_ages)
    }
      
    # output detrended data
    sub_dtrend <- data.frame(site_id = unique(subdat$site_id), site_name = unique(subdat$site_name), entity_id = unique(subdat$entity_id), entity_name = unique(subdat$entity_name), longitude = unique(subdat$longitude), latitude = unique(subdat$latitude),
                             sample_id = subdat$sample_id,
                             interp_age = subdat$interp_age,
                             d18O_detrended = subdat$detrended_d18O)
    
    dtrend_dat <- rbind(dtrend_dat, sub_dtrend)  
  }
}

# no bp's example fig (for Fig. 5)
no_bp_eg <- dtrend_dat %>% filter(entity_id == 590)
p1 <- ggplot(data = no_bp_eg, aes(x = interp_age, y = d18O_detrended)) + geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = seq(7800,8400,100)) +
  scale_y_reverse(breaks = seq(-0.6,0.6,0.2), labels = round(seq(-0.6,0.6,0.2), digits = 1)) +
  xlab("Age (years BP)") + ylab(expression(paste("detrended ", delta^18, "O ", "(\u2030)"))) +
  annotate(geom = "text", x = 7900, y = -0.7, label = "a) no breakpoints (entity id = 590)")

## summarise results
x <- bp_dat %>% group_by(entity_id) %>% summarise(n_bp = n())

## entities with 2 breakpoints
x2 <- x %>% filter(n_bp == 2) # 22

# plot each record individually
for (i in seq(9,27,9)){
  entities <- x2$entity_id[(i-8):i]
  subdat <- dtrend_dat %>% filter(entity_id %in% entities)
  sub_bp <- bp_dat %>% filter(entity_id %in% entities)

  file_no <- ifelse(i == 9, 1, ifelse(i == 18, 2, 3))
  filename <- paste("signals_82_", file_no, ".pdf", sep = "")

  pdf(filename, width = 20/2.54, height = 18/2.54)
  ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
    geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
    geom_vline(data = sub_bp, aes(xintercept = CI2_5), col = "red", lty = 2) +
    geom_vline(data = sub_bp, aes(xintercept = CI97_5), col = "red", lty = 2) +
    facet_wrap(.~ entity_id)
  dev.off()
}

# calculate anomaly and duration for each record
anom_2bp <- data.frame()
for (i in unique(x2$entity_id)){
  sub_bp <- bp_dat %>% filter(entity_id == i)
  subdat <- dtrend_dat %>% filter(entity_id == i)
  
  # sig diff from base:
  subdat <- subdat %>% arrange(interp_age)
  sub_bp <- sub_bp %>% arrange(bp)
  subdat$grp <- cut(subdat$interp_age, breaks = c(min(subdat$interp_age), sub_bp$bp, max(subdat$interp_age)), labels = 1:3)
  t_test <- t.test(x = subdat[which(subdat$grp == 2),"d18O_detrended"], y = subdat[which(subdat$grp %in% c(1,3)),"d18O_detrended"])
  
  # overall
  d18O_event <- subdat %>% filter(interp_age >= min(sub_bp$bp) & interp_age <= max(sub_bp$bp)) %>% summarise(mean(d18O_detrended))
  d18O_base <- subdat %>% filter(interp_age <= min(sub_bp$bp) | interp_age >= max(sub_bp$bp)) %>% summarise(mean(d18O_detrended))
  d18Osd_base <- subdat %>% filter(interp_age <= min(sub_bp$bp) | interp_age >= max(sub_bp$bp)) %>% summarise(sd(d18O_detrended))
  
  #lower limit uncert
  d18O_event_lower <- subdat %>% filter(interp_age >= min(sub_bp$CI2_5) & interp_age <= max(sub_bp$CI97_5)) %>% summarise(mean(d18O_detrended))
  d18O_base_lower <- subdat %>% filter(interp_age <= min(sub_bp$CI2_5) | interp_age >= max(sub_bp$CI97_5)) %>% summarise(mean(d18O_detrended))
  
  # upper limit uncert
  d18O_event_upper <- subdat %>% filter(interp_age >= min(sub_bp$CI97_5) & interp_age <= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
  d18O_base_upper <- subdat %>% filter(interp_age <= min(sub_bp$CI97_5) | interp_age >= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
  
  anom <- d18O_event - d18O_base
  anom_lower <- d18O_event_lower - d18O_base_lower
  anom_upper <- d18O_event_upper - d18O_base_upper
  
  sub_df <- data.frame(unique(sub_bp[,c(1:6)]), 
                       min_bp = as.numeric(min(sub_bp$bp)), max_bp = as.numeric(max(sub_bp$bp)),
                       sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-as.numeric(min(sub_bp$bp)))),"sample_id"],
                       sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-as.numeric(max(sub_bp$bp)))),"sample_id"],
                       anom = as.numeric(anom),
                       #anom_lower = as.numeric(anom_lower),
                       #anom_upper = as.numeric(anom_upper),
                       d18O_base = as.numeric(d18O_base),
                       d18O_event = as.numeric(d18O_event),
                       d18Osd_base = as.numeric(d18Osd_base),
                       ttest_Pval = t_test$p.value)
  
  anom_2bp <- rbind(anom_2bp, sub_df)
}

# filter to those with anomalies greater than the base s.d.
anom_2bp <- anom_2bp %>% mutate(diff_from_base_sd = abs(anom) - abs(d18Osd_base)) 
no_signal <- rbind(no_signal, filter(anom_2bp, diff_from_base_sd <= 0)[,c("site_id","site_name","entity_id","entity_name","longitude","latitude")])
anom_2bp <- anom_2bp %>% filter(diff_from_base_sd > 0)


# filter to those that are sig diff from base
no_signal <- rbind(no_signal, filter(anom_2bp, ttest_Pval >= 0.001)[,c("site_id","site_name","entity_id","entity_name","longitude","latitude")])
anom_2bp <- anom_2bp %>% filter(ttest_Pval <= 0.001)


# 2 bp's example fig (for Fig. 5)
bp2_eg <- dtrend_dat %>% filter(entity_id == 199)
bp2_eg_bps <- bp_dat %>% filter(entity_id == 199)
bp2_eg_basemean <- bp2_eg %>% filter(interp_age >= 7800 & interp_age <= min(bp2_eg_bps$bp) | interp_age >= max(bp2_eg_bps$bp)) %>% summarise(mean(d18O_detrended))
bp2_eg_eventmean <- bp2_eg %>% filter(interp_age >= min(bp2_eg_bps$bp) & interp_age <= max(bp2_eg_bps$bp)) %>% summarise(mean(d18O_detrended))

p2 <- ggplot() + geom_line(data = bp2_eg, aes(x = interp_age, y = d18O_detrended)) +
  geom_line(mapping = aes(x = c(7800, min(bp2_eg_bps$bp)), y = rep(bp2_eg_basemean$`mean(d18O_detrended)`, 2)), col = "dark grey", lty = 5) +
  geom_line(mapping = aes(x = c(8400, max(bp2_eg_bps$bp)), y = rep(bp2_eg_basemean$`mean(d18O_detrended)`, 2)), col = "dark grey", lty = 5) +
  geom_line(mapping = aes(x = c(min(bp2_eg_bps$bp), max(bp2_eg_bps$bp)), y = rep(bp2_eg_eventmean$`mean(d18O_detrended)`)), col = "#00BFC4", lty = 5) +
  geom_vline(data = bp2_eg_bps, aes(xintercept = bp)) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(7800,8400,100)) +
  scale_y_reverse(breaks = seq(-1,1,0.2), labels = round(seq(-1,1,0.2), digits = 1)) +
  xlab("Age (years BP)") + ylab(expression(paste("detrended ", delta^18, "O ", "(\u2030)"))) +
  annotate(geom = "text", x = 7900, y = -0.9, label = "a) 2 breakpoints (entity id = 199)")


### For records with >2 breakpoints, identify the significant segments (grps) (= event)

## entities with 3 breakpoints
x3 <- x %>% filter(n_bp == 3) # 20

# visualise entities
#for (i in seq(9,28,9)){
  entities <- x3$entity_id[(i-8):i]
  subdat <- dtrend_dat %>% filter(entity_id %in% entities)
  sub_bp <- bp_dat %>% filter(entity_id %in% entities)

  file_no <- ifelse(i == 9, 1, ifelse(i == 18, 2, 3))
  filename <- paste("signals_82_3bp_", file_no, ".pdf", sep = "")

  p <- ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
    geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
    facet_wrap(.~ entity_id)
#  pdf(filename, width = 20/2.54, height = 18/2.54)
#  print(p)
#  dev.off()
#}

# 
grps <- data.frame()
for (i in unique(x3$entity_id)){ #for each record
  # filter to record#s breakpoints and detrended data
  subdat <- dtrend_dat %>% filter(entity_id == i) %>% arrange(interp_age)
  sub_bp <- bp_dat %>% filter(entity_id == i) %>% arrange(bp)
  
  # segment record based on breakpoints
  subdat$grp <- with(subdat, ifelse(interp_age <= sub_bp$bp[1] | interp_age >= sub_bp$bp[3], 1,
                                    ifelse(interp_age >= sub_bp$bp[1] & interp_age <= sub_bp$bp[2], 2, 3)))
  subdat$grp <- as.factor(subdat$grp)
  
  # Tukey's HSD of segments
  sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
  sub.av <- aov(sub.lm)
  
  sub_hsd <- TukeyHSD(sub.av)$grp
  
  ## Identifying the event
  sub_hsd2 <- as.data.frame(sub_hsd[which(grepl("1", rownames(sub_hsd))),]) #sig diff of segments from base
  sub_hsd3 <- t(as.data.frame(sub_hsd[which(!grepl("1", rownames(sub_hsd))),])) #sig diff of non-base segments from one another
  
  n_sig <- nrow(sub_hsd2[which(sub_hsd2[,4] <= 0.001),]) #how many segments are sig diff from base?
  
  if (n_sig == 0){ #no segments sig diff from base, i.e. no event
    no_signal <- rbind(no_signal, unique(subdat[,c("site_id","site_name", "entity_id", "entity_name", "longitude", "latitude")]))
  } else if (n_sig == 1){ #only 1 segment is sig diff from base
    if (sub_hsd3[,4] < 0.001){ #non-base segements are sig diff from one another
      event_grp <- rownames(sub_hsd2[which(sub_hsd2[,4] <= 0.001),])
      event_grp <- as.numeric(sub("-1", "", event_grp))
      grps <- rbind(grps, data.frame(entity_id = i, group = event_grp))
    } else { #non-base segments are insig from one another - merge
      grps <- rbind(grps, data.frame(entity_id = rep(i,2), group = c(2,3)))
    }
  } else {
    if (sub_hsd3[,4] < 0.001){ #both non-base segments are sig diff from base
      grp_2 <- mean(subdat[which(subdat$grp == 2),"d18O_detrended"]) - mean(subdat[which(subdat$grp == 1),"d18O_detrended"])
      grp_3 <- mean(subdat[which(subdat$grp == 3),"d18O_detrended"]) - mean(subdat[which(subdat$grp == 1),"d18O_detrended"])
      
      if (all(c(grp_2, grp_3) > 0) | all(c(grp_2, grp_3) < 0)){ # both segments have same direction anomalies - merge
        grps <- rbind(grps, data.frame(entity_id = rep(i, 2), group = c(2,3)))
      } else { # segments have diff direction anomalies - pick segment with the biggest anomalies
        base_mean <- mean(subdat[which(subdat$grp == 1),"d18O_detrended"])
        anom_grp2 <- round(grp_2 - base_mean, digits = 1)
        anom_grp3 <- round(grp_3 - base_mean, digits = 1)
        
        if (abs(anom_grp2) > abs(anom_grp3)){ event_grp <- 2} else if (abs(anom_grp2) < abs(anom_grp3)) { 
          event_grp <- 3 } else {
            event_grp <- rownames(sub_hsd2[which(sub_hsd2[,4] == min(sub_hsd2$`p adj`)),])
            event_grp <- as.numeric(sub("-1", "", event_grp))
          }
        grps <- rbind(grps, data.frame(entity_id = i, group = event_grp))
      }
    }
  }
}



## entities with 4 breakpoints
x4 <- x %>% filter(n_bp == 4) # 7

# visualise
subdat <- dtrend_dat %>% filter(entity_id %in% x4$entity_id)
sub_bp <- bp_dat %>% filter(entity_id %in% x4$entity_id)

#filename <- "signals_82_4bp.pdf"

p <- ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
  facet_wrap(.~ entity_id)
#pdf(filename, width = 20/2.54, height = 18/2.54)
#print(p)
#dev.off()

# grps
for (i in unique(x4$entity_id)){
  subdat <- dtrend_dat %>% filter(entity_id == i) %>% arrange(interp_age)
  sub_bp <- bp_dat %>% filter(entity_id == i) %>% arrange(bp)
  
  subdat$grp <- with(subdat, ifelse(interp_age <= sub_bp$bp[1] | interp_age >= sub_bp$bp[4], 1,
                                    ifelse(interp_age >= sub_bp$bp[1] & interp_age <= sub_bp$bp[2], 2, 
                                           ifelse(interp_age >= sub_bp$bp[2] & interp_age <= sub_bp$bp[3], 3, 4))))
  subdat$grp <- as.factor(subdat$grp)
  sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
  sub.av <- aov(sub.lm)
  
  sub_hsd <- TukeyHSD(sub.av)$grp
  #sub_hsd2 <- HSD.test(sub.av, trt = 'grp', alpha = 0.001)
  
  ## QC checks
  sub_hsd2 <- as.data.frame(sub_hsd[which(grepl("1", rownames(sub_hsd))),])
  sub_hsd3 <- as.data.frame(sub_hsd[which(!grepl("1", rownames(sub_hsd))),])
  
  n_sig <- nrow(sub_hsd2[which(sub_hsd2[,4] <= 0.001),])
  
  if (n_sig == 0){ # 220
    no_signal <- rbind(no_signal, unique(subdat[,c("site_id","site_name", "entity_id", "entity_name", "longitude", "latitude")]))
    
  } else if (n_sig == 1){ #244
    event_grp <- rownames(sub_hsd2[which(sub_hsd2[,4] <= 0.001),])
    event_grp <- as.numeric(sub("-1", "", event_grp))
    grps <- rbind(grps, data.frame(entity_id = i, group = event_grp))
    
  } else if (n_sig == 2){
    if (nrow(sub_hsd3[which(sub_hsd3$`p adj` <= 0.001),]) == 0){ # all events are insig from one another, ent = 52
      grps <- rbind(grps, data.frame(entity_id = rep(i, 3), group = 2:4))
    } else { #327, 351, 442
      event_grp <- rownames(sub_hsd2[which(sub_hsd2[,4] <= 0.001),])
      event_grp <- as.numeric(sub("-1","", event_grp))
      grps <- rbind(grps, data.frame(entity_id = rep(i, length(event_grp)+1), group = event_grp[1]:event_grp[2]))
    }
    
  } else if (n_sig == 3){
    if (nrow(sub_hsd3[which(sub_hsd3$`p adj` <= 0.001),]) == 3){
      grp_2 <- mean(subdat[which(subdat$grp == 2),"d18O_detrended"])
      grp_3 <- mean(subdat[which(subdat$grp == 3),"d18O_detrended"])
      grp_4 <- mean(subdat[which(subdat$grp == 4),"d18O_detrended"])
      
      if (all(c(grp_2, grp_3, grp_4) > 0) | all(c(grp_2, grp_3, grp_4) < 0)){ #254
        grps <- rbind(grps, data.frame(entity_id = rep(i, 3), group = 2:4))
      }
    } else {
      event_grp <- rownames(sub_hsd3[which(sub_hsd3[,4] >= 0.001),])
      event_grp <- as.numeric(as.vector(strsplit(event_grp, "-"))[[1]])
      grps <- rbind(grps, data.frame(entity_id = rep(i, length(event_grp)+1), group = event_grp[1]:event_grp[2]))
    }
  }
}


# 3 bp's example fig (for Fig. 5)
bp3_eg <- dtrend_dat %>% filter(entity_id == 129)
bp3_eg_bps <- bp_dat %>% filter(entity_id == 129)
bp3_eg_basemean <- bp3_eg %>% filter(interp_age >= 7800 & interp_age <= min(bp3_eg_bps$bp) | interp_age >= max(bp3_eg_bps$bp)) %>% summarise(mean(d18O_detrended))
bp3_eg_eventmean1 <- bp3_eg %>% filter(interp_age >= min(bp3_eg_bps$bp) & interp_age <= bp3_eg_bps$bp[2]) %>% summarise(mean(d18O_detrended))
bp3_eg_eventmean2 <- bp3_eg %>% filter(interp_age >= bp3_eg_bps$bp[2] & interp_age <= bp3_eg_bps$bp[3]) %>% summarise(mean(d18O_detrended))

bp3_eg_bps <- bp3_eg_bps %>% arrange(bp)

p3 <- ggplot() + geom_line(data = bp3_eg, aes(x = interp_age, y = d18O_detrended)) +
  geom_line(mapping = aes(x = c(7800, min(bp3_eg_bps$bp)), y = rep(bp3_eg_basemean$`mean(d18O_detrended)`, 2)), col = "dark grey", lty = 5) +
  geom_line(mapping = aes(x = c(8400, max(bp3_eg_bps$bp)), y = rep(bp3_eg_basemean$`mean(d18O_detrended)`, 2)), col = "dark grey", lty = 5) +
  geom_line(mapping = aes(x = c(min(bp3_eg_bps$bp), bp3_eg_bps$bp[2]), y = rep(bp3_eg_eventmean1$`mean(d18O_detrended)`)), col = "#00BFC4", lty = 5) +
  geom_line(mapping = aes(x = c(bp3_eg_bps$bp[2], bp3_eg_bps$bp[3]), y = rep(bp3_eg_eventmean2$`mean(d18O_detrended)`,2)), col = "#00BFC4", lty = 5) +
  geom_vline(data = bp3_eg_bps, aes(xintercept = bp)) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(7800,8400,100)) +
  scale_y_reverse(breaks = seq(-1.4,1,0.2), labels = round(seq(-1.4,1,0.2), digits = 1)) +
  xlab("Age (years BP)") + ylab(expression(paste("detrended ", delta^18, "O ", "(\u2030)"))) +
  annotate(geom = "text", x = 7900, y = -1.3, label = "a) >2 breakpoints (entity id = 129)")



## 5 bp's
x5 <- x %>% filter(n_bp == 5) # 3

# visualise
subdat <- dtrend_dat %>% filter(entity_id %in% x5$entity_id)
sub_bp <- bp_dat %>% filter(entity_id %in% x5$entity_id)

#filename <- "signals_82_5bp.pdf"

p <- ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
  facet_wrap(.~ entity_id)
#pdf(filename, width = 20/2.54, height = 18/2.54)
#print(p)
#dev.off()

subdat <- dtrend_dat %>% filter(entity_id == 591)
sub_bp <- bp_dat %>% filter(entity_id == 591)

#subdat$grp <- with(subdat, ifelse(interp_age <= sub_bp$bp[1] | interp_age >= sub_bp$bp[5], 1,
#                                  ifelse(interp_age >= sub_bp$bp[1] & interp_age <= sub_bp$bp[2], 2, 
#                                         ifelse(interp_age >= sub_bp$bp[2] & interp_age <= sub_bp$bp[3], 3,
#                                                ifelse(interp_age >= sub_bp$bp[3] & interp_age <= sub_bp$bp[4], 4, 5)))))
#subdat$grp <- as.factor(subdat$grp)
subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp, max(interp_age)), labels = 1:6))

sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
sub.av <- aov(sub.lm)

sub_hsd <- TukeyHSD(sub.av)$grp

grps <- rbind(grps, data.frame(entity_id = rep(591,2), group = c(3,4)))


# >2 bp's, no signal, example fig (for Fig. 5)
bp5_eg <- dtrend_dat %>% filter(entity_id == 305)
bp5_eg_bps <- bp_dat %>% filter(entity_id == 305)
bp5_eg_basemean <- bp5_eg %>% filter(interp_age >= 7800 & interp_age <= min(bp5_eg_bps$bp) | interp_age >= max(bp5_eg_bps$bp)) %>% summarise(mean(d18O_detrended))
bp5_eg_eventmean1 <- bp5_eg %>% filter(interp_age >= min(bp5_eg_bps$bp) & interp_age <= bp5_eg_bps$bp[2]) %>% summarise(mean(d18O_detrended))
bp5_eg_eventmean2 <- bp5_eg %>% filter(interp_age >= bp5_eg_bps$bp[2] & interp_age <= bp5_eg_bps$bp[3]) %>% summarise(mean(d18O_detrended))
bp5_eg_eventmean3 <- bp5_eg %>% filter(interp_age >= bp5_eg_bps$bp[3] & interp_age <= bp5_eg_bps$bp[4]) %>% summarise(mean(d18O_detrended))
bp5_eg_eventmean4 <- bp5_eg %>% filter(interp_age >= bp5_eg_bps$bp[4] & interp_age <= bp5_eg_bps$bp[5]) %>% summarise(mean(d18O_detrended))


bp5_eg_bps <- bp5_eg_bps %>% arrange(bp)

p4 <- ggplot() + geom_line(data = bp5_eg, aes(x = interp_age, y = d18O_detrended)) +
  geom_line(mapping = aes(x = c(7800, min(bp5_eg_bps$bp)), y = rep(bp5_eg_basemean$`mean(d18O_detrended)`, 2)), col = "dark grey", lty = 5) +
  geom_line(mapping = aes(x = c(8400, max(bp5_eg_bps$bp)), y = rep(bp5_eg_basemean$`mean(d18O_detrended)`, 2)), col = "dark grey", lty = 5) +
  geom_line(mapping = aes(x = c(min(bp5_eg_bps$bp), bp5_eg_bps$bp[2]), y = rep(bp5_eg_eventmean1$`mean(d18O_detrended)`)), col = "#00BFC4", lty = 5) +
  geom_line(mapping = aes(x = c(bp5_eg_bps$bp[2], bp5_eg_bps$bp[3]), y = rep(bp5_eg_eventmean2$`mean(d18O_detrended)`)), col = "#00BFC4", lty = 5) +
  geom_line(mapping = aes(x = c(bp5_eg_bps$bp[3], bp5_eg_bps$bp[4]), y = rep(bp5_eg_eventmean3$`mean(d18O_detrended)`)), col = "#00BFC4", lty = 5) +
  geom_line(mapping = aes(x = c(bp5_eg_bps$bp[4], bp5_eg_bps$bp[5]), y = rep(bp5_eg_eventmean4$`mean(d18O_detrended)`)), col = "#00BFC4", lty = 5) +
  geom_vline(data = bp5_eg_bps, aes(xintercept = bp)) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(7800,8400,100)) +
  scale_y_reverse(breaks = seq(-0.6,1,0.2), labels = round(seq(-0.6,1,0.2), digits = 1)) +
  xlab("Age (years BP)") + ylab(expression(paste("detrended ", delta^18, "O ", "(\u2030)"))) +
  annotate(geom = "text", x = 7900, y = -0.5, label = "a) >2 breakpoints (entity id = 305)")

## Fig. 5
library(cowplot)
pdf("Fig5.pdf", width = 18/2.54, height = 12/2.54)
plot_grid(p1,p2,p3,p4, ncol = 2, align = "hv")
dev.off()


# calculate anomalies
anom_bp <- data.frame()
x_all <- rbind(x3,x4,x5)
for (i in unique(grps$entity_id)){
  sub_bp <- bp_dat %>% filter(entity_id == i) %>% arrange(bp)
  #sub_bp <- sub_bp[order(sub_bp$bp),]
  subdat <- dtrend_dat %>% filter(interp_age >= 7400 & interp_age <= 9000 & 
                                    entity_id == i)
  
  if (nrow(sub_bp) == 3){ lab_vec = 1:4 } else if (nrow(sub_bp) == 4) { lab_vec = 1:5 } else { lab_vec = 1:6}
  subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp,max(interp_age)), labels = lab_vec))
  
  event_grp <- as.numeric(grps[which(grps$entity_id == i),2])
  event_grp <- na.omit(event_grp)
  
  event_d18O <- subdat %>% filter(grp %in% event_grp) %>% summarise(mean_d18O = mean(d18O_detrended))
  base_d18O <- subdat %>% filter(grp %in% c(min(lab_vec), max(lab_vec))) %>% summarise(mean_d18O = mean(d18O_detrended))
  d18Osd_base <- subdat %>% filter(grp %in% c(min(lab_vec), max(lab_vec))) %>% summarise(d18Osd_base = sd(d18O_detrended))
  #base_d18O <- subdat %>% filter(interp_age >= 7400 & interp_age <= 7900 | interp_age >= 8500 & interp_age <= 9000) %>% summarise(mean_d18O = mean(d18O_detrended))
  x_sub <- data.frame(unique(subdat[which(subdat$entity_id == i),1:6]),
                      min_bp = sub_bp$bp[min(event_grp,na.rm = T)-1], max_bp = sub_bp$bp[max(event_grp, na.rm = T)], 
                      sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[min(event_grp,na.rm = T)-1])),"sample_id"],
                      sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[max(event_grp, na.rm = T)])),"sample_id"],
                      anom = as.numeric(event_d18O) - as.numeric(base_d18O),
                      d18O_base = as.numeric(base_d18O), d18O_event = as.numeric(event_d18O),
                      d18Osd_base = as.numeric(d18Osd_base$d18Osd_base))
  anom_bp <- rbind(anom_bp, x_sub)
}

# filter to those with anomalies greater than the base s.d.
anom_bp <- anom_bp %>% mutate(diff_from_base_sd = abs(anom) - abs(d18Osd_base)) 
no_signal <- rbind(no_signal, filter(anom_bp, diff_from_base_sd < 0)[,c("site_id","site_name","entity_id","entity_name","longitude","latitude")])
anom_bp <- anom_bp %>% filter(diff_from_base_sd >= 0)


# combine
all_dat <- rbind(anom_2bp[,-c(14:17)], anom_bp[,c(1:13)])

# calculate duration
all_dat$duration <- all_dat$max_bp - all_dat$min_bp
all_dat$event_centre <- all_dat$max_bp - (all_dat$duration/2)

## remove spurious results: where there are >1 records from the same site with disagreement, select the records with higher resolution
unique_sites <- all_dat %>% group_by(site_id, site_name, latitude, longitude) %>% summarise(n())
mult_sites <- unique_sites %>% filter(`n()` >= 2)

mult_out <- data.frame()
for (i in unique(mult_sites$site_id)){
  if (i == 8){ next }
  subdat <- all_dat %>% filter(site_id == i)
  subres <- res_out %>% filter(site_id == i)
  
  if (nrow(subdat) == 2){
    if (all(subdat$anom >= 0) | all(subdat$anom <= 0)){
      mult_out <- rbind(mult_out, subdat)
    } else {
      highres_ent <- subres %>% filter(mean_res == min(subres$mean_res))
      mult_out <- rbind(mult_out, filter(subdat, entity_id == highres_ent$entity_id))
    }
  } else {
    if (all(subdat$anom >= 0) | all(subdat$anom <= 0)){
      mult_out <- rbind(mult_out, subdat)
    } else {
      lowres_ent <- subres %>% filter(mean_res == max(subres$mean_res))
      mult_out <- rbind(mult_out, filter(subdat, entity_id != lowres_ent$entity_id))
    }
  }
}

all_dat <- all_dat %>% filter(!site_id %in% mult_out$site_id)
all_dat <- rbind(all_dat, mult_out)

## add non-SISAL records
nonSISAL <- read.csv("C:/Users/ph805612/OneDrive - University of Reading/Documents/abrupt_Holocene/nonSISAL_82_signals.csv")

all_dat <- rbind(all_dat, nonSISAL)


## save
write.csv(all_dat, "spel_82_signals.csv", row.names = F)
write.csv(no_signal, "C:/Users/ph805612/OneDrive - University of Reading//Documents/abrupt_Holocene/spel_nosignal_82.csv", row.names = F)


### Make table of results for supplement
supp_out <- all_dat %>% select(entity_id, site_name, longitude, latitude, max_bp, min_bp, duration, anom) %>% 
  arrange(entity_id, site_name)
colnames(supp_out)[5:6] <- c("start","end")

write.csv(supp_out, "anom_table.csv", row.names = F)


# significant difference between regions?

Europe <- all_dat %>% filter(latitude >= 20 & longitude >= -10 & longitude <= 45)
Asia <- all_dat %>% filter(latitude >= 0 & latitude <= 45 & longitude >= 50 & longitude <= 150)
S_America <- all_dat %>% filter(latitude >= -30 & latitude <= 0 & longitude >= -100 & longitude <= -30)
global <- all_dat 

Europe$region <- "Europe"
Asia$region <- "Asia"
S_America$region <- "S_America"
global$region <- "global"

all_dat2 <- rbind(Europe, Asia, global)


t.test(x = Asia$max_bp, y = Europe$max_bp)
t.test(x = Asia$min_bp, y = Europe$min_bp)
t.test(x = abs(Europe$anom), y = abs(Asia$anom))
