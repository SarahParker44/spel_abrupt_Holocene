## Constrain 10.2-10.7 kyr signals

setwd("C:/Users/sarah/OneDrive - University of Reading/Documents/abrupt_Holocene/speleothem_8_2_kyr_signals/")
setwd("C:/Users/ph805612/OneDrive - University of Reading/Documents/abrupt_Holocene/")

library(dplyr)
library(strucchange)
library(ggplot2)
library(RMySQL)

# connect to SISAL database
mydb <- dbConnect(MySQL(), user = "root", password = "BevRed921", dbname = "sisal_v2", 
                  host = "localhost")

# load data
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING (site_id) JOIN sample USING (entity_id) JOIN original_chronology USING (sample_id) JOIN d18O USING (sample_id)
                       WHERE (interp_age BETWEEN 10100 AND 10800);")
Raw_Data <- Raw_Data %>% filter(entity_status != "superseded")

write.csv(Raw_Data, "rawdat_102_107.csv", row.names = F)
