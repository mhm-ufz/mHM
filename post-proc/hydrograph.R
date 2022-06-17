######################################################################################################### 
##                            ---------------------------------------------------------------------------
## ========================== Hydrograph Generation from mHM streamflow output (discharge.nc)
##                            ----------------------------------------------------------------------------
## ---------- Code developer: 
## -------------------------  Pallav Kumar Shrestha pallav-kumar.shrestha@ufz.de ------
## -------------------------  07 October 2021 ----------------------------------------
#########################################################################################################

#### Gives one set of skill scores for the whole simulation period

# (1.1) LIBRARIES PREPARATION ===================

# Check for the required packages
list.of.packages <- c("ggplot2", "hydroGOF", "ncdf4", "chron", "xts", "reshape", 
                      "stringr", "ragg")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# Open libraries/ packages
# silence on
oldw <- getOption("warn")
options(warn = -1) 

# Open libraries/ packages
library(ggplot2)
library(hydroGOF) # for using fucntions KGE and NSE
library(ncdf4) 
library(chron)
library(xts) 
library(reshape) # for melt
library(stringr) # for str_pad
library(ragg)    # for scaling plot

options(warn = oldw) 
# silence off



plot_hydrograph <- function(path, gID, title_text, fType, ylimit ){


  # ========  CONTROL  =============================================
  
  # Parameters
  misVal = -9999.0
  
  # Graph control
  use_labels <- c("observed", "mHM")
  use_colors <- c("red", "blue")
  use_linetypes <- c(2, 1)
  
  
  # ========  READ  =============================================
  
  # Source file
  if (missing(fType)){
    fName = "discharge.nc"
  } else {
    if (fType == "txt"){
      fName = "daily_discharge.out"
    } else {
      fName = "discharge.nc"
    }
  }

  if (fName == "discharge.nc"){

    # Read the netCDF discharge file
    ncin <- nc_open(paste(path,fName,sep = "/"))
    # get VARIABLES
    q_obs <- ncvar_get(ncin, paste("Qobs_", str_pad(gID, 10, pad = "0"), sep = ""))
    q_sim <- ncvar_get(ncin, paste("Qsim_", str_pad(gID, 10, pad = "0"), sep = ""))
    # Read time attribute
    nctime <- ncvar_get(ncin,"time")
    tunits <- ncatt_get(ncin,"time","units")
    # Close file
    nc_close(ncin)
  
    # Prepare the time origin
    tustr <- unlist(strsplit(tunits$value, " "))
    tdstr <- unlist(strsplit((rev(tustr))[2], "-"))
    tmonth <- as.integer(tdstr)[2]
    tday <- as.integer(tdstr)[3]
    tyear <- as.integer(tdstr)[1]
    tchron <- chron(dates. = nctime/24, origin=c(tmonth, tday, tyear)) # nctime (hours)

  } else {

    # Reading the TEXT discharge file
    data = data.frame(read.delim(paste(path,fName,sep="/"), header = TRUE, sep = ""))  # reading all the data
    q_obs <- data[[paste("Qobs_", str_pad(gID, 10, pad = "0"), sep = "")]]
    q_sim <- data[[paste("Qsim_", str_pad(gID, 10, pad = "0"), sep = "")]]
    nData <- length(data[,1])
    dStart <- as.Date(paste(data[1,4],"-",data[1,3],"-",data[1,2],sep=""))  # Infering the start date
    dEnd <- as.Date(paste(data[nData,4],"-",data[nData,3],"-",data[nData,2],sep=""))  # Infering the end date
    tchron <- seq.Date(dStart,dEnd, by= "days")

  }
  
  
  
  
  
  # ========  PROCESS  =============================================
  
  # Replacing missing values by NA
  q_obs[q_obs == misVal] <- NA
  q_sim[q_sim == misVal] <- NA
  
  # convert to xts
  q_obs <- xts(as.numeric(q_obs), order.by = tchron) # xts/ time series object created
  q_sim <- xts(as.numeric(q_sim), order.by = tchron)
  
  
  # Bind
  q <- cbind(q_obs, q_sim)
  # Melt
  q_df <- data.frame(q)
  q_df$id <- rownames(q_df)
  q_melted <- melt(q_df, measure.vars=c("q_obs", "q_sim"))
  # id must be date class
  q_melted$id <- rep(seq.Date(as.Date(tchron[1]),as.Date(tail(tchron, n = 1)), by= "days"), 2)
  
  
  # Preparing METRICS annotations for the graph
  statKGE <- round(KGE(q_sim,q_obs,na.rm = TRUE),2)  # KGE 
  statNSE <- round(NSeff(q_sim,q_obs,na.rm = TRUE),2)  # NSE
  statPBS <- round(pbias(q_sim,q_obs,na.rm = TRUE),2)  # PBIAS

  statPosX1 <- 
             # start date
             as.Date(tchron[1]) +
             # date window
             ( as.Date(tail(tchron, n = 1)) - as.Date(tchron[1]) ) *
             # start date as fraction of window
             0.09
  statPosX2 <- 
            # start date
            as.Date(tchron[1]) +
            # date window
            ( as.Date(tail(tchron, n = 1)) - as.Date(tchron[1]) ) *
            # start date as fraction of window
            0.25
  qmax <- max(q_obs, q_sim, na.rm = TRUE)
  
  # y limit
  if (missing(ylimit)){
    mylimit = NULL
  } else {
    mylimit = c(0, ylimit)
    qmax = ylimit
  }
  
  statPosY1 <- qmax*0.9   # determining position for statistics
  statPosY2 <- qmax*0.8
  statPosY3 <- qmax*0.7
  statPosY4 <- qmax*0.6
  
  translucentPosY1 <- qmax*0.55
  translucentPosY2 <- qmax*0.95
  translucentPosX1 <- as.Date(tchron[1])
  translucentPosX2 <- as.Date(tail(tchron, n = 1))
  
  # adjust pdf width for long hydrographs
  nyrs = as.numeric(format(tail(tchron, n = 1),'%Y')) - as.numeric(format(tchron[1],'%Y')) + 1
  if (nyrs > 10){
    mydatebreaks = "5 years"
    pdf_width = 12 + (nyrs - 10)
  } else {
    mydatebreaks = "1 year"
    pdf_width = 12
  }
  
  
  # ========  PLOT  =============================================
  
  date_range = paste(as.numeric(format(tchron[1],'%Y')), as.numeric(format(tail(tchron, n = 1),'%Y')), sep = "-")
  caption_text <- paste("KGE: ", statKGE, "     NSE: ", statNSE, "     PBIAS: ", statPBS, " %")
  
  # Plotting the hydrograph
  
  main <- ggplot() +
    # hydrographs
    geom_line(data = q_melted, aes( x = id, y = value, color = as.factor(variable), linetype = as.factor(variable) ), size = 1, alpha = 1) +
    
    scale_color_manual(values = use_colors,
                       labels = use_labels) +
    
    scale_linetype_manual(values = use_linetypes,
                          labels = use_labels) +
    
    labs(title = title_text, caption = caption_text) +
    
    theme(
      text=element_text(family = "Helvetica", colour = "black"),
      axis.ticks.length=unit(-0.2, "cm"),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.text.x = element_text(size=14, margin = margin(t = 10), colour = "black"),
      axis.title.x = element_text(size=16, margin = margin(t = 10), colour = "black"),
      axis.text.y = element_text(size=14, margin = margin(r = 10), colour = "black"),
      axis.title.y.left  = element_text(size=16, margin = margin(r = 15), colour = "black", hjust = c(0.5)),
      axis.title.y.right = element_blank(),
      plot.title = element_text(size = 16, colour = "black", hjust = c(0), margin = margin(b = -50, t = 10), face = "bold"),
      plot.caption = element_text(size = 16, colour = "black", hjust = c(1)),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.background = element_blank(),
      panel.grid.major = element_line(colour = alpha("black", 0.5), size=0.2, linetype = 3),
      legend.position = "top",
      legend.key = element_blank(),
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(1.5, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.box.margin = margin(t = 20, unit = "pt"),
      legend.text = element_text(size=16, colour = "black", hjust = c(0), margin = margin(r = 30, unit = "pt")),
      legend.title = element_blank(),
      legend.background = element_blank()) +
    
    scale_x_date(name = "Time", date_breaks= mydatebreaks, date_labels = "%Y", expand = c(0,0)) + # duplicating the axis for the top was not possible with date axis
  
    scale_y_continuous(name = expression(paste("Streamflow [",m^{3},".",s^{-1},"]")), limits = mylimit , sec.axis = dup_axis(name ="", labels = c()), expand = c(0,0))  # adding extra space at the top for annotations
  
  # Output
  ggsave(main, file=paste(path,"/",gID,"_hydrograph.png",sep=""), width = 18, height = 5, units = "in")

}
