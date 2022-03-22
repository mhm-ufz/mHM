###################################################################################################
##                   -------------------------------------------
## ==================   Plot Script for mHM's discharge file    ===================================
##                   -------------------------------------------
##   
##  Author:    Pallav Kumar Shrestha    (pallav-kumar.shrestha@ufz.de)
##             07 October 2021
##
##  Usage:     Place this file alongside discharge.nc file and type the following in
##             the command line -
##
##                              Rscript hydrograph_call.R
##
##  Output:    <stationid>_hydrograph.png
##
##  Detail:    Creates hydrograph plot of the gauge specified in this file
##
##  Reference: 
##
##  Modifications:
##
###################################################################################################


# Source taylor-made functions
source("hydrograph.R")


## Variables
path <- "." # path to the discharge.nc file. 
station_name <- "The Interesting Gauge" # Name of the station/ gauge
station_id <- "398" # Station ID (should be same as in mhm.nml, do not worry about the padded zeroes)

## Generate texts
title_text <- paste("Gauge: ", station_name, " . ", station_id) # Plot title 
  
## Call the hydrograph function 
plot_hydrograph( path, station_id, title_text, "txt" ) # txt - read from daily_dicharge.txt, nc - read from discharge.nc

