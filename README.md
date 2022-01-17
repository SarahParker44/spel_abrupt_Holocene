# spel_abrupt_Holocene

This github repository provides the code used to produce results and figures from Parker, S.E. and Harrison, S.P. The timing, duration and magnitude of the 8.2 ka event in global speleothem records, *in progress*

Contact: Sarah Parker, s.parker@pgr.reading.ac.uk

## Data used
The code in this repository uses the SISAL (Speleothem Isotopes Synthesis and Analysis) database, available from XX. Database is in MySQL and R connects to this using the RMySQL package. 

## Structure of the respository
### Holocene abrupt event detection
"Holocene_abrupt_events.R" - Detection of abrupt event through the Holocene
"abrupt_event_significance.R" - Determines the signficance of Holocene bins using randomly generated time series with red noise. Plots Fig. 2.

### 8.2 ka event analysis
"signals_8_2_kyr.R" - determines presence and timing of an 8.2 ka isotope signal in records, then calculates the duration and anomalies for these excursions. Calculates the global, European and Asian  median timing, duration and absolute magnitude. Calculaes whether there is a significant difference between regions using t-tests. Plots Fig. 5. 
"plot_maps.R" - plots anomalies on map using results from "signals_8_2_kyr.R". Also plots data from Table. S3 and temperature and precipitation anomalies from Morrill et al. 2013 (Proxy benchmarks for intercomparison of 8.2 ka simulations. Clim. Past 9, 423–432). Plots Fig. 4. 

### Figure of records using in this study
"plot_fig1_map.R" - plots records using the Holoene abrupt event analysis and 8.2 ka event analysis. Plots Fig. 1. 

### Analysis of other Holocene periods of abrupt climate change
"supplement/..." - files calculate the abrupt events in significant Holocene periods. Plots the timing and anomalies for each period. Used to plot Fig. S1. 
