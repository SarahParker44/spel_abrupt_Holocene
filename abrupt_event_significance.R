setwd("C:/Users/sarah/OneDrive - University of Reading/Documents/abrupt_Holocene/")

library(dplyr)
library(strucchange)
library(ggplot2)

## load dtrend_dat
dtrend_dat <- read.csv("Hol_dtrend_dat.csv")

# remove overlapping data
df <- data.frame(x = seq(300,12000,300),
                  y = c(0,0,500,500,1000,1000,1500,1500,2000,2000,2500,3000,3000,3500,3500,4000,4500,4500,5000,5000,5500,6000,6000,6500,7000,7000,7500,7500,8000,8000,8500,9000,9000,9500,9500,10000,10500,10500,11000,11000))
df_no_rep <- data.frame()
for (i in unique(dtrend_dat$entity_id)){
  subdat <- dtrend_dat %>% filter(entity_id == i)
  
  #select ages with >1 row
  rep_ages <- subdat %>% group_by(interp_age) %>% summarise(n()) %>% filter(`n()` > 1)
  
  # 
  x <- subdat %>% filter(interp_age %in% rep_ages$interp_age)
  #
  subdat <- subdat %>% filter(!interp_age %in% rep_ages$interp_age)
  
  y <- data.frame()
  for (j in rep_ages$interp_age){
    xx <- x %>% filter(interp_age == j)
    
    xx$bin <- 300 * ceiling(xx$interp_age/300)
    #select one
    xx <- xx[which(xx$win_start == (df[which(df$x == unique(xx$bin)),]$y)),]
    
    y = rbind(y, xx)
  }
  subdat <- rbind(subdat, y[,-11])
  subdat <- subdat %>% arrange(interp_age)
  
  df_no_rep <- rbind(df_no_rep, subdat)
}

# 300 year bins
df_no_rep$bin <- cut(df_no_rep$interp_age, breaks = seq(0,12000,300), labels = seq(300,12000,300))

df_no_rep <- na.omit(df_no_rep)

entities_exclude <- df_no_rep %>% group_by(entity_id, bin) %>% summarise(n()) %>% filter(`n()` <10)

df_no_rep <- df_no_rep %>% group_by(entity_id, bin) %>% mutate(n()) %>% filter(`n()` >10)

## mean npoints per bin
npoints <- df_no_rep %>% group_by(entity_id, bin) %>% summarise(n_points = n()) %>% group_by(bin) %>% summarise(mean_n_points = mean(n_points))

# mean ar (autoregression coefficient) per bin
get_ar_coef <- function(x){
  xx <- arima(x, c(1,0,0))
  
  return(xx$coef[1])
}
ar_coef <- df_no_rep %>% group_by(entity_id, bin) %>% summarise(ar_coef = get_ar_coef(d18O_detrended))
ar_coef <- ar_coef %>% group_by(bin) %>% summarise(mean_ar = mean(ar_coef))


# for each bin randomaly generate red-noise time series with the same npoints and ar as the actual data (1000 times)
bin_pcent <- data.frame()
for (i in seq(300,12000,300)){
  bin_npoints <- npoints %>% filter(bin == i)
  bin_ar <- ar_coef %>% filter(bin == i)
  
  set.seed(123)
  mat <- replicate(1000, arima.sim(n = bin_npoints$mean_n_points, model = list(order = c(1,0,0), ar = bin_ar$mean_ar)))
  
  out_vec <- c()
  for (j in 1:1000){
    bp <- breakpoints(mat[,j] ~ 1)
    out_vec[j] <- length(bp$breakpoints)
  }
  
  bin_pcent <- rbind(bin_pcent, data.frame(bin = i, pcent = length(out_vec[which(out_vec >= 2)])/1000))
}

# reshape df
xx <- bin_pcent

bin_pcenta <- bin_pcent; bin_pcentb <- bin_pcent
bin_pcenta$bin_age <- bin_pcenta$bin-300
bin_pcentb$bin_age <- bin_pcentb$bin
bin_pcent <- rbind(bin_pcenta, bin_pcentb)

bin_pcent <- bin_pcent %>% arrange(bin_age, bin)

bp_pcent <- read.csv("Hol_bp_pcent.csv")

# plot fig. 2
pdf("Hol_plot.pdf", width = 13/2.54, height = 8/2.54)
ggplot() + 
  geom_line(data = bp_pcent, aes(x = bin_age, y = pcent*100)) + 
  geom_line(data = bin_pcent, aes(x = bin_age, y = pcent*100), col = "red") +
  
  geom_segment(aes(x = 8200, y = 80, xend = 8200, yend = 73), arrow = arrow(length = unit(0.1, "cm")), col = "red") +
  geom_text(aes(x = 8200, y = 82, label = "8.2 ka"), col = "red") +
  scale_x_continuous(breaks = seq(0,12000,1000), expand = c(0.01,0.01)) +
  ylab("% entities") + xlab("Age (years BP)") +
  theme_bw()
dev.off()

