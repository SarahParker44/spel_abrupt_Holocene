## Constrain 1.5-1.8 kyr signals

setwd(".../spel_abrupt_Holocene")

library(dplyr)
library(strucchange)
library(ggplot2)
library(RMySQL)

# connect to SISAL database
mydb <- dbConnect(MySQL(), user = "root", password = "BevRed921", dbname = "sisal_v2", 
                  host = "localhost")

# load 1.4 to 1.9 ka  data
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING (site_id) JOIN sample USING (entity_id) JOIN original_chronology USING (sample_id) JOIN d18O USING (sample_id)
                       WHERE (interp_age BETWEEN 1400 AND 1900);")
Raw_Data <- Raw_Data %>% filter(entity_status != "superseded")

#write.csv(Raw_Data, "rawdat_15_18.csv", row.names = F)

## Select entities with a sufficient resolution between 7800 and 8400 years (of 30 years)
min_res <- 30

# load function for calculating mean temporal resolution (excluding hiatuses and gaps) - 'get_ent_sampling'
source("speleothem_8_2_kyr_signals/entity_sampling_mean_res.R")

res_out <- data.frame()
for (i in unique(Raw_Data$entity_id)){ # for every entity
  subdat <- Raw_Data %>% filter(entity_id == i & interp_age >= 1400 & interp_age <= 1900)
  if (nrow(subdat) <= 1){ next }
  mean_res <- get_ent_sampling(entity_id = i, age_start = 1400, age_end = 1900)$sampling_mean
  
  res_out <- rbind(res_out,
                   data.frame(unique(subdat[c("site_id","site_name","entity_id")]),
                              mean_res = mean_res))
}

res_out <- res_out %>% filter(mean_res <= 30)

## filter to records of a sufficient resolution
res_out$record_length <- NA
for (i in 1:length(res_out$entity_id)){
  ent_id <- res_out[i,"entity_id"]
  subdat <- Raw_Data %>% filter(entity_id == ent_id & interp_age >=1400 & interp_age <= 1900)
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
  subdat <- Raw_Data %>% filter(entity_id == i & interp_age >= 1400 & interp_age <= 1900)
  
  if (nrow(subdat) <= 15){ next }
  
  ## detrend
  subdat_lm <- lm(d18O_measurement ~ interp_age, data = subdat)
  lm_predicted <- predict(subdat_lm)
  subdat$detrended_d18O <- residuals(subdat_lm)
  
  bp <- breakpoints(subdat$detrended_d18O ~ 1) #breakpoint analysis
  #bpts_sum <- summary(bp) 
  #opt_brks <- opt_bpts(bpts_sum$RSS["BIC",]) #optimal no. breakpoints
  
  #if (length(opt_brks) > 1){ 
  #  opt_brks <- as.numeric(names(which.min(bpts_sum$RSS["BIC",]))) } 
  #else 
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
      
      # output detrended data
      sub_dtrend <- data.frame(site_id = unique(subdat$site_id), site_name = unique(subdat$site_name), entity_id = unique(subdat$entity_id), entity_name = unique(subdat$entity_name), longitude = unique(subdat$longitude), latitude = unique(subdat$latitude),
                               sample_id = subdat$sample_id,
                               interp_age = subdat$interp_age,
                               d18O_detrended = subdat$detrended_d18O)
      
      dtrend_dat <- rbind(dtrend_dat, sub_dtrend)
    }
  }
}

## summarise results
x <- bp_dat %>% group_by(entity_id) %>% summarise(n_bp = n())

## entities with 2 breakpoints
x2 <- x %>% filter(n_bp == 2) # 33

i=36
entities <- x2$entity_id[(i-8):i]
subdat <- dtrend_dat %>% filter(entity_id %in% entities)
sub_bp <- bp_dat %>% filter(entity_id %in% entities)
ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
  geom_vline(data = sub_bp, aes(xintercept = CI2_5), col = "red", lty = 2) +
  geom_vline(data = sub_bp, aes(xintercept = CI97_5), col = "red", lty = 2) +
  facet_wrap(.~ entity_id)

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
anom_2bp <- anom_2bp %>% filter(ttest_Pval <= 0.001) #21

ggplot(data = anom_2bp, aes(x = longitude, y = latitude, fill = anom)) +
  geom_point(shape = 21, size = 3) + 
  borders("world") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue")

## entities with 3 breakpoints
x3 <- x %>% filter(n_bp == 3) # 12

# visualise entities
#for (i in seq(9,28,9)){
i=18
entities <- x3$entity_id[(i-8):i]
subdat <- dtrend_dat %>% filter(entity_id %in% entities)
sub_bp <- bp_dat %>% filter(entity_id %in% entities)
ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
  facet_wrap(.~ entity_id)

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
  
  if (n_sig == 0){ #no segments sig diff from base, i.e. no event ## 538
    no_signal <- rbind(no_signal, unique(subdat[,c("site_id","site_name", "entity_id", "entity_name", "longitude", "latitude")]))
    print(paste("no signal =", i))
  } else if (n_sig == 1){ #only 1 segment is sig diff from base
    if (sub_hsd3[,4] < 0.001){ #non-base segements are sig diff from one another
      event_grp <- rownames(sub_hsd2[which(sub_hsd2[,4] <= 0.001),])
      event_grp <- as.numeric(sub("-1", "", event_grp))
      grps <- rbind(grps, data.frame(entity_id = i, group = event_grp))
      print(paste("grp 1 =",i))
    } else { #non-base segments are insig from one another - merge
      grps <- rbind(grps, data.frame(entity_id = rep(i,2), group = c(2,3)))
      print(paste("grp 2 =", i))
    }
  } else {
    if (sub_hsd3[,4] < 0.001){ #both non-base segments are sig diff from base
      grp_2 <- mean(subdat[which(subdat$grp == 2),"d18O_detrended"]) - mean(subdat[which(subdat$grp == 1),"d18O_detrended"])
      grp_3 <- mean(subdat[which(subdat$grp == 3),"d18O_detrended"]) - mean(subdat[which(subdat$grp == 1),"d18O_detrended"])
      
      if (all(c(grp_2, grp_3) > 0) | all(c(grp_2, grp_3) < 0)){ # both segments have same direction anomalies - merge
        grps <- rbind(grps, data.frame(entity_id = rep(i, 2), group = c(2,3)))
        print(paste("grp 3 =", i))
      } else { # segments have diff direction anomalies - pick segment with the biggest anomalies
        base_mean <- mean(subdat[which(subdat$grp == 1),"d18O_detrended"])
        anom_grp2 <- round(grp_2 - base_mean, digits = 1)
        anom_grp3 <- round(grp_3 - base_mean, digits = 1)
        
        if (abs(anom_grp2) > abs(anom_grp3)){ 
          event_grp <- 2
        } else if (abs(anom_grp2) < abs(anom_grp3)) { 
          event_grp <- 3 
        } else {
          event_grp <- rownames(sub_hsd2[which(sub_hsd2[,4] == min(sub_hsd2$`p adj`)),])
          event_grp <- as.numeric(sub("-1", "", event_grp))
        }
        grps <- rbind(grps, data.frame(entity_id = i, group = event_grp))
        print(paste("grp 4 = ",i))
      }
    } else { # non-base segments are insig from one another
      grps <- rbind(grps, data.frame(entity_id = rep(i,2), group = 2:3))
    }
  }
}

## entities with 4 breakpoints
x4 <- x %>% filter(n_bp == 4) # 6

# visualise
subdat <- dtrend_dat %>% filter(entity_id %in% x4$entity_id)
sub_bp <- bp_dat %>% filter(entity_id %in% x4$entity_id)
ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
  facet_wrap(.~ entity_id)

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
  
  if (n_sig == 0){ 
    no_signal <- rbind(no_signal, unique(subdat[,c("site_id","site_name", "entity_id", "entity_name", "longitude", "latitude")]))
    print(paste("no signal = ",i))
  } else if (n_sig == 1){ 
    event_grp <- rownames(sub_hsd2[which(sub_hsd2[,4] <= 0.001),])
    event_grp <- as.numeric(sub("-1", "", event_grp))
    grps <- rbind(grps, data.frame(entity_id = i, group = event_grp))
    print(paste("grp1 =", i))
  } else if (n_sig == 2){
    if (nrow(sub_hsd3[which(sub_hsd3$`p adj` <= 0.001),]) == 0){ # all events are insig from one another
      grps <- rbind(grps, data.frame(entity_id = rep(i, 3), group = 2:4))
    } else {
      event_grp <- rownames(sub_hsd2[which(sub_hsd2[,4] <= 0.001),])
      event_grp <- as.numeric(sub("-1","", event_grp))
      grps <- rbind(grps, data.frame(entity_id = rep(i, length(event_grp[1]:event_grp[2])), group = event_grp[1]:event_grp[2]))
      print(paste("grp 2 = ", i))
    }
    
  } else if (n_sig == 3){
    if (nrow(sub_hsd3[which(sub_hsd3$`p adj` <= 0.001),]) == 3){
      grp_2 <- mean(subdat[which(subdat$grp == 2),"d18O_detrended"])
      grp_3 <- mean(subdat[which(subdat$grp == 3),"d18O_detrended"])
      grp_4 <- mean(subdat[which(subdat$grp == 4),"d18O_detrended"])
      
      if (all(c(grp_2, grp_3, grp_4) > 0) | all(c(grp_2, grp_3, grp_4) < 0)){ 
        grps <- rbind(grps, data.frame(entity_id = rep(i, 3), group = 2:4))
        print(paste("grp 3 =",i))
      } else {
        no_signal <- rbind(no_signal, unique(subdat[,c("site_id","site_name", "entity_id", "entity_name", "longitude", "latitude")]))
        print(paste("no signal = ",i))
      }
    } else {
      event_grp <- rownames(sub_hsd3[which(sub_hsd3[,4] >= 0.001),])
      event_grp <- as.numeric(as.vector(strsplit(event_grp, "-"))[[1]])
      grps <- rbind(grps, data.frame(entity_id = rep(i, length(event_grp)+1), group = event_grp[1]:event_grp[2]))
    }
  }
}

# calculate anomalies
anom_bp <- data.frame()
x_all <- rbind(x3,x4)
for (i in unique(grps$entity_id)){
  sub_bp <- bp_dat %>% filter(entity_id == i) %>% arrange(bp)
  #sub_bp <- sub_bp[order(sub_bp$bp),]
  subdat <- dtrend_dat %>% filter(interp_age >= 1400 & interp_age <= 1900 & 
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
anom_bp <- anom_bp %>% filter(diff_from_base_sd >= 0)  ##18->10


# combine
all_dat <- rbind(anom_2bp[,-c(14:17)], anom_bp[,c(1:13)])

# calculate duration
all_dat$duration <- all_dat$max_bp - all_dat$min_bp
all_dat$event_centre <- all_dat$max_bp - (all_dat$duration/2)

# remove spurious results: where there are >1 records from the same site with disagreement, select the records with higher resolution
unique_sites <- all_dat %>% group_by(site_id, site_name, latitude, longitude) %>% summarise(n())
mult_sites <- unique_sites %>% filter(`n()` >= 2)

mult_out <- data.frame()
for (i in unique(mult_sites$site_id)){
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

##              
ggplot(data = all_dat, aes(x = longitude, y = latitude, fill = anom)) +
  geom_point(shape = 21, size = 3) + 
  borders("world") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue")  

## little regional coherence

# sites with signal verus sites with no signal
nrow(all_dat %>% group_by(site_id) %>% summarise(n())) #27
nrow(no_signal %>% group_by(site_id) %>% summarise(n())) #42

27/(27+42) #39% (just above threshold)

## save dat
write.csv(all_dat, "anom_1500_1900.csv", row.names = F)
write.csv(no_signal, "nosignal_1500_1900.csv", row.names = F)

