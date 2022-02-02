###################################################################################################
##                   -------------------------------------------
## ==================   Animator for mHM output netCDF files    ===================================
##                   -------------------------------------------
##   
##  Author:    Pallav Kumar Shrestha    (pallav-kumar.shrestha@ufz.de)
##             10 May 2019
##
##  Usage:     Place this file alongside *_Fluxes_States.nc file/s and type the following in
##             the command line -
##
##                              Rscript animate1.R
##
##  Output:    animate.gif, animate.pdf
##
##  Detail:    Creates multi-plots of all variables in mHM and mRM flux states files
##
##  Reference: https://sites.google.com/view/climate-access-cooperative/code?authuser=0 
##
##  Modifications:
##             Pallav Shrestha - 27 May 2019
##             "the plot size and font size now scale with number of variables being plot"
##
##             Pallav Shrestha - 2 June 2019
##             "NCL colors introduced for better visualization"
##
###################################################################################################

cat("
============================================================
NOTE:
      At the moment, animate1.R behaves well when mHM netCDF 
      output files have lesser time slices. E.g. MONTHLY or
      yearly output or daily output containing few months
============================================================\n\n")


#============================================================================= (1) PREPARATION ====

# (1.1) LIBRARIES PREPARATION ===================

# Check for the required packages
list.of.packages <- c("ggplot2", "ncdf4", "graphics", "RColorBrewer", "chron", "plyr", 
                      "lattice", "animation", "gridExtra", "grid", "classInt")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# Open libraries/ packages
# silence on
oldw <- getOption("warn")
options(warn = -1) 

library(ggplot2)
library(ncdf4) 
library(graphics) 
library(RColorBrewer)
library(chron)
library(plyr)
library(lattice)
library(animation)
library(gridExtra)  # for using "grid.arrange" function
library(grid)
library(classInt)

options(warn = oldw) 
# silence off

# (1.2) FILE READ ===============================

# Input file names
fName = c("mRM_Fluxes_States.nc", "mHM_Fluxes_States.nc" )
nfiles = 2

# Check and Open the netCDF files
nvar=0
if (!file.exists(fName[1])) {
  nfiles <- nfiles - 1
  print("Warning: mRM output file is missing")
} else {
  ncin_mrm <- nc_open(fName[1])
  ncin_common <- ncin_mrm
  nvar <- nvar + 1
}
if (!file.exists(fName[2])) {
  nfiles <- nfiles - 1
  print("Warning: mHM output file is missing")
} else {
  ncin_mhm <- nc_open(fName[2])
  ncin_common <- ncin_mhm
  varlist <- names(ncin_mhm$var)
  varlist <- varlist[ varlist != "time_bnds" & varlist !="lat" & varlist != "lon"] # excluding time and lat lon variables
  nvar <- nvar + length(varlist)
}
if (nfiles == 0) {
  print("Error : both output files are missing. Nothing to plot!")
  q() # quit the program!
}



# Read time attribute
nctime <- ncvar_get(ncin_common,"time")
tunits <- ncatt_get(ncin_common,"time","units") 
nt <- dim(nctime)


# (1.3) COLORS & PLOT PREPARATION ========================

# Prepare the time origin
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth <- as.integer(unlist(tdstr)[2])
tday <- as.integer(unlist(tdstr)[3])
tyear <- as.integer(unlist(tdstr)[1])
tfinal <- chron(nctime/24, origin=c(tmonth, tday, tyear)) # nctime (hours) is converted to days

# GIF speed setting based on time resolution
temp_res_check <- tfinal[2] - tfinal[1]
if (temp_res_check <= 1) {
  gifinterval <- 0.2  # daily
} else if (temp_res_check >= 32) {
  gifinterval <- 2    # yearly
} else {
  gifinterval <- 1    # monthly
}

# Determine the matrix layout of multiplot
cols <- ceiling(sqrt(nvar))
rows <- floor(sqrt(nvar)) +
              as.integer( nvar/ ( ceiling(sqrt(nvar))*floor(sqrt(nvar)) + 1) )

# Lat-Lon suitable intervals for snapping
latlon_cutpts <- c(0.1, 0.2, 0.25, 0.5, 1, 2, 5)

# COLOR PALETTES from NCL stored as COLOR MATRIX
#     (source: https://www.ncl.ucar.edu/Document/Graphics/ColorTables/Files/)

# silence on
oldw <- getOption("warn")
options(warn = -1) 
clr_mat   <-    # 1, cmocean_tempo, 256
                c("#FFF6F4","#FDF5F3","#FCF4F1","#FBF3F0","#F9F2EE","#F8F1ED","#F7F0EB","#F5EFEA","#F4EEE8","#F2EDE7","#F1ECE5","#F0EBE4","#EEEAE2","#EDEAE1","#EBE9DF","#EAE8DE","#E9E7DD","#E7E6DB","#E6E5DA","#E4E4D8","#E3E3D7","#E2E2D6","#E0E2D4","#DFE1D3","#DDE0D1","#DCDFD0","#DBDECF","#D9DDCD","#D8DDCC","#D6DCCB","#D5DBC9","#D3DAC8","#D2D9C7","#D1D8C5","#CFD8C4","#CED7C3","#CCD6C1","#CBD5C0","#C9D4BF","#C8D4BE","#C6D3BC","#C5D2BB","#C3D1BA","#C2D1B9","#C0D0B7","#BFCFB6","#BDCEB5","#BCCEB4","#BACDB3","#B9CCB2","#B7CBB0","#B6CBAF","#B4CAAE","#B3C9AD","#B1C8AC","#B0C8AB","#AEC7AA","#ACC6A9","#ABC5A8","#A9C5A6","#A8C4A5","#A6C3A4","#A4C3A3","#A3C2A2","#A1C1A1","#A0C0A0","#9EC09F","#9CBF9F","#9BBE9E","#99BE9D","#97BD9C","#96BC9B","#94BC9A","#92BB99","#91BA98","#8FBA97","#8DB997","#8BB896","#8AB795","#88B794","#86B693","#85B593","#83B592","#81B491","#7FB390","#7DB390","#7CB28F","#7AB18E","#78B18E","#76B08D","#74AF8D","#72AF8C","#71AE8B","#6FAD8B","#6DAD8A","#6BAC8A","#69AB89","#67AB89","#65AA88","#63A988","#61A987","#5FA887","#5DA786","#5BA686","#59A685","#57A585","#56A485","#54A484","#52A384","#50A284","#4EA183","#4BA183","#49A083","#479F82","#459F82","#439E82","#419D82","#3F9C81","#3D9C81","#3B9B81","#3A9A81","#389981","#369880","#349880","#329780","#309680","#2E9580","#2C947F","#2A937F","#29937F","#27927F","#25917F","#24907F","#228F7E","#218E7E","#1F8D7E","#1E8D7E","#1C8C7E","#1B8B7D","#1A8A7D","#19897D","#17887D","#16877C","#16867C","#15857C","#14847C","#13847B","#13837B","#12827B","#12817B","#11807A","#117F7A","#117E7A","#117D79","#117C79","#117B79","#117A78","#117978","#117878","#117777","#117677","#127676","#127576","#127476","#137375","#137275","#137174","#147074","#146F73","#146E73","#156D73","#156C72","#166B72","#166A71","#166971","#176870","#176770","#17666F","#18656F","#18656E","#18646E","#19636D","#19626D","#19616C","#19606C","#1A5F6B","#1A5E6B","#1A5D6A","#1A5C6A","#1A5B69","#1B5A68","#1B5968","#1B5867","#1B5867","#1B5766","#1B5666","#1C5565","#1C5465","#1C5364","#1C5263","#1C5163","#1C5062","#1C4F62","#1C4E61","#1C4D61","#1C4C60","#1C4C5F","#1C4B5F","#1C4A5E","#1C495E","#1C485D","#1C475D","#1C465C","#1C455B","#1C445B","#1C435A","#1C425A","#1C4259","#1C4158","#1C4058","#1B3F57","#1B3E57","#1B3D56","#1B3C56","#1B3B55","#1B3A54","#1B3954","#1B3853","#1A3753","#1A3652","#1A3651","#1A3551","#1A3450","#1A3350","#19324F","#19314F","#19304E","#192F4D","#192E4D","#182D4C","#182C4C","#182B4B","#182A4B","#18294A","#17284A","#172749","#172648","#172548","#172447","#162347","#162246","#162146","#162045","#151F45","#151E44","#151D44")
clr_mat   <-    # 2, MPL_PuBu, 128
                rbind(clr_mat, c("#FFF7FC","#FEF6FB","#FDF5FB","#FCF4FA","#FBF3F9","#F9F2F9","#F8F1F8","#F7F0F8","#F6EFF7","#F5EFF7","#F3EDF6","#F3EDF6","#F1EBF5","#F0EBF5","#EFE9F4","#EEE9F3","#ECE7F3","#EAE6F2","#E9E5F1","#E7E3F0","#E5E2F0","#E3E0EF","#E2E0EE","#E0DDED","#DEDCEC","#DCDBEC","#DBDAEB","#D9D8EA","#D7D6E9","#D5D5E9","#D4D4E8","#D1D2E7","#CFD1E6","#CDD0E6","#CACFE5","#C7CDE4","#C5CCE4","#C3CBE3","#BFC9E2","#BDC8E2","#BAC7E1","#B7C6E0","#B5C4E0","#B2C3DF","#AFC2DE","#AEC1DE","#AABFDD","#A7BEDC","#A4BDDB","#A1BCDB","#9EBADA","#9BB9D9","#98B8D8","#96B7D8","#92B5D7","#8EB4D6","#8BB3D5","#88B2D5","#85B0D4","#82AFD3","#7FAED2","#7DADD2","#78ABD1","#75AAD0","#72A8CF","#6EA7CE","#6AA5CD","#66A4CC","#62A2CB","#5EA1CA","#5A9FC9","#569DC8","#529CC8","#4E9AC7","#4A99C6","#4998C5","#4396C4","#3F94C3","#3B92C2","#3791C1","#348FC0","#318DBF","#2E8BBE","#2A89BD","#2787BC","#2485BB","#2183BA","#1E81B9","#1B7FB8","#187DB7","#157BB6","#137AB5","#0F77B4","#0C75B3","#0873B2","#0571B1","#056FAF","#056EAD","#056CAA","#056BA8","#056AA6","#0568A4","#0567A2","#05669F","#04649D","#04639B","#046199","#046198","#045F94","#045D92","#045C90","#045A8E","#04588A","#045687","#045484","#045280","#03507D","#034E7A","#034B76","#034973","#034770","#03456C","#034369","#034267","#023F62","#023D5F","#023A5C","#023858"))
clr_mat   <-    # 3, NEO_div_vegetation_a, 256
                rbind(clr_mat, c("#530000","#540000","#540000","#550000","#560000","#560000","#570000","#580000","#5A0000","#5A0000","#5B0000","#5B0000","#5D0000","#5E0000","#5F0000","#600000","#610000","#610000","#620000","#630000","#640000","#650000","#670000","#670000","#690000","#6A0000","#6B0000","#6C0000","#6D0000","#6E0000","#700000","#710000","#720000","#740000","#740000","#760000","#770000","#790000","#790000","#7B0000","#7C0000","#7D0000","#7F0000","#800000","#830401","#830702","#860B04","#870D05","#8A1106","#8C1507","#8E1908","#901D0A","#921F0B","#93220C","#96270D","#972A0E","#992E0F","#9C3211","#9E3512","#9F3913","#A13C14","#A23F15","#A54417","#A74818","#A94B19","#AB501A","#AD531B","#B0571D","#B25B1E","#B35D1F","#B66321","#B86722","#BB6C24","#BD6F25","#BF7326","#C0762A","#C1782D","#C2792F","#C37C33","#C57F37","#C6813A","#C7843E","#C88641","#C98945","#CB8C49","#CC8E4C","#CC904F","#CE9252","#CF9656","#D09859","#D19B5D","#D29D60","#D49F64","#D5A368","#D6A56B","#D7A86F","#D9AA73","#D9AC75","#DBB07A","#DCB27D","#DDB581","#DFB885","#E0BA88","#E1BD8D","#E2C090","#E3C293","#E5C597","#E6C79A","#E7CA9E","#E9CDA3","#EAD0A6","#EBD2A9","#ECD5AC","#EDD7B0","#EFDAB4","#F0DDB8","#F1E1BD","#F2E3C0","#F3E5C3","#F5E8C7","#F6EBCB","#F7EDCE","#F8EFD1","#FAF4D6","#FBF5D9","#FCF8DD","#FDFBE1","#FFFFE5","#FEFFE5","#FAFDE0","#F8FBDD","#F6FBDB","#F4F9D7","#F1F8D4","#EFF7D1","#ECF6CD","#EAF4CA","#E7F3C6","#E5F2C3","#E2F1BF","#E0EFBB","#DDEEB9","#DAEDB5","#D8ECB1","#D6EAAF","#D2E9AB","#D0E8A8","#CDE7A5","#CBE5A1","#C8E49F","#C6E39B","#C3E297","#C1E194","#BFDF91","#BCDE8D","#BADD8B","#B7DB86","#B5DA83","#B3D981","#B0D87D","#AED77A","#ABD576","#A9D574","#A6D371","#A3D26D","#A1D16A","#9ED066","#9CCE64","#9ACD60","#97CC5D","#96CB5B","#93CA57","#90C954","#8EC750","#8CC74E","#89C54A","#87C448","#85C344","#83C242","#80C03D","#7FC03C","#7BBE39","#77BC37","#72B935","#6EB733","#6BB532","#67B331","#63B02F","#5FAE2D","#5BAC2C","#58AB2B","#54A829","#51A628","#4DA426","#4AA325","#46A124","#429F22","#409D20","#3B9A1E","#38991D","#34971C","#31951A","#2D9319","#299118","#278F17","#238D15","#1F8B13","#1B8913","#198711","#168510","#13840F","#0F810D","#0C800C","#0C7E0C","#0C7D0C","#0B7B0B","#0B790B","#0B780B","#0B770B","#0B750B","#0B750B","#0A730A","#0A710A","#0A700A","#0A6E0A","#0A6D0A","#0A6C0A","#0A6B0A","#096909","#096709","#096609","#096609","#096409","#096309","#086108","#086108","#085F08","#085F08","#085D08","#085C08","#085A08","#085A08","#075907","#075807","#075607","#075507","#075507","#075407","#075207","#075107","#065106","#065006","#064F06","#064E06","#064D06","#064C06"))
clr_mat   <-    # 4, GMT_drywet, 60
                rbind(clr_mat, c("#8C672D","#967133","#A17B39","#AB853F","#B68F44","#C09A4A","#CAA450","#D5AE56","#DFB85C","#EAC361","#ECCA66","#E6CE6A","#E0D26D","#DBD571","#D5D974","#CFDD78","#C9E17B","#C3E57F","#BDE982","#B8ED86","#AEEF8D","#A1EF97","#94EFA1","#87EFAB","#7AEFB5","#6DEFBF","#60EFC9","#53EFD3","#46EFDD","#39EFE7","#30E9EC","#2CDDEC","#29D1ED","#25C5ED","#21BAED","#1DAEEE","#19A2EE","#1696EE","#128AEE","#0E7EEF","#0D72EC","#1067E7","#135BE1","#154FDC","#1843D6","#1A37D1","#1D2BCB","#201FC6","#2213C0","#2507BA","#2504B4","#2209AD","#1F0EA6","#1C139F","#191898","#161D91","#13228A","#102783","#0D2C7C","#0A3175"))
clr_mat   <-    # 5, NCV_blu_red, 256
                rbind(clr_mat, c("#2A007F","#2A0082","#2A0085","#2A0088","#2A008B","#29008E","#290091","#290094","#280097","#28009A","#28009D","#2700A0","#2700A3","#2600A6","#2500A9","#2500AC","#2400AF","#2300B2","#2200B5","#2200B8","#2100BB","#2000BE","#1F00C1","#1E00C4","#1C00C7","#1B00CA","#1A00CD","#1900D0","#1800D3","#1600D6","#1500D9","#1300DC","#1200DF","#1000E2","#0F00E5","#0D00E8","#0C00EB","#0A00EE","#0800F1","#0600F4","#0400F7","#0200FA","#0000FD","#0102FF","#0407FF","#070CFF","#0A11FF","#0D16FF","#101AFF","#131FFF","#1624FF","#1928FF","#1C2DFF","#1F31FF","#2236FF","#253AFF","#283FFF","#2B43FF","#2E47FF","#314CFF","#3450FF","#3754FF","#3A58FF","#3D5CFF","#4060FF","#4364FF","#4668FF","#496CFF","#4C70FF","#4F73FF","#5277FF","#557BFF","#587FFF","#5B82FF","#5E86FF","#6189FF","#648DFF","#6790FF","#6A93FF","#6D97FF","#709AFF","#739DFF","#76A0FF","#79A3FF","#7CA6FF","#7FAAFF","#82ACFF","#85AFFF","#88B2FF","#8BB5FF","#8EB8FF","#91BBFF","#94BDFF","#97C0FF","#9AC3FF","#9DC5FF","#A0C8FF","#A3CAFF","#A6CDFF","#A9CFFF","#ACD1FF","#AFD3FF","#B2D6FF","#B5D8FF","#B8DAFF","#BBDCFF","#BEDEFF","#C1E0FF","#C4E2FF","#C7E4FF","#CAE6FF","#CDE8FF","#D0E9FF","#D3EBFF","#D6EDFF","#D9EEFF","#DCF0FF","#DFF1FF","#E2F3FF","#E5F4FF","#E8F6FF","#EBF7FF","#EEF8FF","#F1FAFF","#F4FBFF","#F7FCFF","#FAFDFF","#FDFEFF","#FFFEFD","#FFFDFA","#FFFCF7","#FFFBF4","#FFFAF1","#FFF8EE","#FFF7EB","#FFF6E8","#FFF4E5","#FFF3E2","#FFF1DF","#FFF0DC","#FFEED9","#FFEDD6","#FFEBD3","#FFE9D0","#FFE8CD","#FFE6CA","#FFE4C7","#FFE2C4","#FFE0C1","#FFDEBE","#FFDCBB","#FFDAB8","#FFD8B5","#FFD6B2","#FFD3AF","#FFD1AC","#FFCFA9","#FFCDA6","#FFCAA3","#FFC8A0","#FFC59D","#FFC39A","#FFC097","#FFBD94","#FFBB91","#FFB88E","#FFB58B","#FFB288","#FFAF85","#FFAC82","#FFA97F","#FFA67C","#FFA379","#FFA076","#FF9D73","#FF9A70","#FF976D","#FF936A","#FF9067","#FF8D64","#FF8961","#FF865E","#FF825B","#FF7F58","#FF7B55","#FF7752","#FF734F","#FF704C","#FF6C49","#FF6846","#FF6443","#FF6040","#FF5C3D","#FF583A","#FF5437","#FF5034","#FF4C31","#FF472E","#FF432B","#FF3F28","#FF3A25","#FF3622","#FF311F","#FF2D1C","#FF2819","#FF2416","#FF1F13","#FF1A10","#FF160D","#FF110A","#FF0C07","#FF0704","#FF0201","#FD0000","#FA0002","#F70004","#F40006","#F10008","#EE000A","#EB000C","#E8000D","#E5000F","#E20010","#DF0012","#DC0013","#D90015","#D60016","#D30018","#D00019","#CD001A","#CA001B","#C7001C","#C4001E","#C1001F","#BE0020","#BB0021","#B80022","#B50022","#B20023","#AF0024","#AC0025","#A90025","#A60026","#A30027","#A00027","#9D0028","#9A0028","#970028","#940029","#910029","#8E0029","#8B002A","#88002A","#85002A","#82002A","#7F002A"))
clr_mat   <-    # 6, MPL_RdBu, 128
                rbind(clr_mat, c("#6A0120","#700321","#760521","#7C0722","#820923","#880A24","#8E0C25","#940E26","#9A1027","#9D1128","#A51429","#A8152A","#B1182B","#B31A2C","#B82230","#B92531","#BD2D35","#C03338","#C23639","#C63E3D","#C94440","#CC4A43","#CD4D44","#D15548","#D45B4B","#D7604D","#D8634F","#DC6B56","#DE715A","#E0765E","#E17960","#E58166","#E7866B","#EA8C6F","#EC9173","#EE9777","#F19C7B","#F29F7D","#F5A784","#F6AB89","#F7AF8F","#F7B394","#F8B89A","#F9BC9F","#F9C0A4","#FAC2A7","#FBC9AF","#FCCDB5","#FCD1BA","#FDD5C0","#FEDAC5","#FEDDCA","#FDDFCD","#FDE0CF","#FCE4D5","#FCE6D9","#FBE8DD","#FBEAE0","#FAECE4","#FAEFE8","#FAF1EC","#F9F2EE","#F9F5F3","#F8F7F7","#F6F7F8","#F3F5F7","#F0F4F6","#EDF3F6","#EAF1F5","#E7F0F5","#E4EEF4","#E1EDF4","#DEECF3","#DBEAF3","#D8E9F2","#D6E8F2","#D2E6F1","#CDE3F0","#C8E1EE","#C3DEED","#BEDCEB","#B9D9EA","#B4D7E8","#AFD4E7","#AAD2E6","#A5CFE4","#A0CDE3","#9BCAE1","#96C8E0","#91C5DE","#8BC1DC","#88BFDB","#7EB9D8","#78B5D6","#72B1D4","#6CADD2","#65A9CF","#5FA5CD","#59A1CB","#539DC9","#4D99C7","#4696C5","#4292C3","#3F8EC1","#3D8BBF","#3A87BD","#3784BC","#3682BB","#327DB8","#2F79B6","#2D75B4","#2A72B3","#276EB1","#246BAF","#2267AD","#1F63A8","#1D5FA2","#1B5B9C","#195697","#175291","#144E8B","#124A85","#10457F","#0F437C","#0C3D73","#09396D","#073467","#053061"))
clr_mat   <-    # 7, MPL_RdYlBu, 128
                rbind(clr_mat, c("#A70126","#AB0526","#AF0926","#B30D26","#B71126","#BB1426","#BF1826","#C31C26","#C72026","#C92226","#CE2726","#D02927","#D62F27","#D83127","#DB382B","#DC3A2C","#E0422F","#E24731","#E34932","#E75036","#E95538","#EB5A3A","#EC5C3B","#F0633E","#F26841","#F46D43","#F56F44","#F67747","#F77C4A","#F7814C","#F8844D","#F98C51","#F99153","#FA9656","#FB9B58","#FCA05A","#FCA55D","#FDA85E","#FEAF62","#FEB366","#FEB769","#FEBB6D","#FEBF71","#FEC374","#FEC778","#FEC97A","#FECF7F","#FED383","#FED787","#FEDA8B","#FEDE8E","#FFE292","#FFE496","#FFE597","#FFE99D","#FFEBA1","#FFEEA4","#FFF0A8","#FFF3AC","#FFF5B0","#FFF8B3","#FFF9B5","#FFFCBB","#FFFFBE","#FEFFC3","#FBFEC7","#F9FDCC","#F6FCD0","#F4FBD5","#F1FAD9","#EFF9DE","#EDF8E2","#EAF7E7","#E8F6EB","#E5F5EF","#E4F5F2","#E0F3F8","#DCF1F7","#D8EFF6","#D4EDF5","#D0EBF4","#CCE9F3","#C7E7F1","#C3E5F0","#BFE3EF","#BBE1EE","#B7DFED","#B2DDEB","#AEDBEA","#AAD8E9","#A6D5E7","#A4D3E6","#9DCEE3","#99CBE1","#94C7DF","#90C4DE","#8CC0DC","#87BDDA","#83B9D8","#7FB6D6","#7AB2D4","#76AFD2","#72ABD0","#6EA7CE","#6BA2CC","#679EC9","#6399C7","#6197C6","#5C91C2","#588CC0","#5488BE","#5183BC","#4D7FB9","#497AB7","#4676B5","#4471B2","#426CB0","#4067AD","#3F62AB","#3D5DA9","#3C58A6","#3A53A4","#394FA1","#384CA0","#35459C","#34409A","#323B98","#313695"))
clr_mat   <-    # 8, MPL_BrBG, 128
                rbind(clr_mat, c("#573105","#5B3406","#5F3706","#643906","#683C07","#6D3E07","#714108","#754408","#7A4608","#7C4809","#834B09","#854D09","#8B510A","#8E520B","#94580F","#965A10","#9C5F14","#A06317","#A26519","#A86B1D","#AC6F1F","#B07222","#B27424","#B87A28","#BC7E2A","#C0822D","#C18430","#C58C3A","#C79140","#CA9646","#CB994A","#CFA053","#D1A559","#D4AA60","#D6B066","#D9B56C","#DBBA72","#DDBC76","#E0C47F","#E2C784","#E4C98A","#E6CC8F","#E8CF95","#E9D29A","#EBD5A0","#ECD7A3","#EFDBAB","#F1DEB0","#F2E1B6","#F4E4BB","#F6E7C1","#F7E9C6","#F7EACA","#F7EBCC","#F7ECD2","#F7EED5","#F7EFD9","#F6F0DD","#F6F1E1","#F6F2E5","#F6F3E9","#F6F3EB","#F6F5F1","#F6F6F5","#F3F5F5","#F0F4F4","#ECF4F2","#E8F3F1","#E5F2F0","#E1F1EF","#DEF0ED","#DAEFEC","#D6EEEB","#D3EEEA","#CFEDE8","#CDECE8","#C8EBE6","#C2E9E3","#BDE6E0","#B7E4DD","#B1E2DB","#ACE0D8","#A6DDD5","#A1DBD2","#9BD9CF","#95D6CC","#90D4CA","#8AD2C7","#85D0C4","#7FCDC1","#79C8BD","#76C6BB","#6DC0B5","#67BCB1","#61B7AD","#5CB3A9","#56AFA5","#50ABA1","#4AA69D","#44A299","#3E9E95","#389A92","#33968E","#2F928A","#2B8E86","#278A82","#23867E","#21847C","#1B7F76","#167B73","#12776F","#0E736B","#0A6F67","#066B63","#02675F","#01645C","#016158","#015D54","#015A51","#01574D","#01534A","#005046","#004D42","#004B40","#00463B","#004337","#004034","#003C30"))
clr_mat   <-    # 9, precip_diff_12lev, 13
                rbind(clr_mat, c("#B66A28","#CD853F","#E1A564","#F5CD84","#F5E09E","#FFF5BA","#FFFFFF","#CDFFCD","#99F0B2","#53BD9F","#6EAAC8","#0570B0","#023858"))
clr_mat   <-    # 10, MPL_YlGnBu, 128
                rbind(clr_mat, c("#FFFFD8","#FEFFD6","#FDFED3","#FCFED1","#FAFECE","#F9FDCB","#F8FDC9","#F7FCC6","#F6FCC4","#F5FCC3","#F4FBBF","#F3FBBE","#F1FABA","#F1FAB9","#EFF9B5","#EEF9B4","#ECF8B1","#EAF7B1","#E8F7B2","#E5F5B2","#E3F4B2","#E0F3B2","#DFF3B2","#DBF1B3","#D9F0B3","#D7EFB3","#D5EFB3","#D2EEB3","#CFEDB4","#CDECB4","#CCEBB4","#C8EAB4","#C4E8B4","#C0E7B5","#BBE5B5","#B7E3B6","#B2E1B6","#B0E0B6","#A9DEB7","#A5DCB8","#A0DAB8","#9CD8B8","#97D7B9","#93D5B9","#8ED3BA","#8CD2BA","#85D0BB","#80CEBB","#7CCCBC","#78CBBC","#75C9BD","#71C8BD","#6DC7BE","#6BC6BE","#65C4BF","#61C2C0","#5DC1C0","#59BFC1","#55BEC1","#51BCC2","#4DBBC2","#4BBAC3","#46B8C4","#42B7C4","#3FB4C4","#3DB2C4","#3BB0C4","#38ADC3","#36ABC3","#34A9C3","#31A6C3","#2FA4C2","#2DA2C2","#2B9FC2","#289DC2","#279CC1","#2498C1","#2296C1","#1F94C1","#1D92C0","#1D8EBF","#1D8BBE","#1D88BC","#1E85BA","#1E82B9","#1E7EB7","#1F7BB6","#1F78B4","#1F75B3","#2072B1","#206EB0","#206DAF","#2168AD","#2165AB","#2161AA","#225EA8","#225CA7","#2259A6","#2256A5","#2254A3","#2351A2","#234EA1","#234C9F","#23499E","#23469D","#23449C","#24419A","#24409A","#243C98","#243997","#243795","#253494","#233291","#21318D","#1F2F89","#1D2E85","#1C2D81","#1A2B7E","#182A7A","#162876","#142772","#12256F","#11246B","#102369","#0D2163","#0B205F","#091E5C","#081D58"))
clr_mat   <-    # 11, NCV_jet, 256
                rbind(clr_mat, c("#00007F","#000083","#000087","#00008B","#00008F","#000093","#000097","#00009B","#00009F","#0000A3","#0000A7","#0000AB","#0000AF","#0000B3","#0000B7","#0000BB","#0000BF","#0000C3","#0000C7","#0000CB","#0000CF","#0000D3","#0000D7","#0000DB","#0000DF","#0000E3","#0000E7","#0000EB","#0000EF","#0000F3","#0000F7","#0000FB","#0000FF","#0004FF","#0008FF","#000CFF","#0010FF","#0014FF","#0018FF","#001CFF","#0020FF","#0024FF","#0028FF","#002CFF","#0030FF","#0034FF","#0038FF","#003CFF","#0040FF","#0044FF","#0048FF","#004CFF","#0050FF","#0054FF","#0058FF","#005CFF","#0060FF","#0064FF","#0068FF","#006CFF","#0070FF","#0074FF","#0078FF","#007CFF","#0080FF","#0084FF","#0088FF","#008CFF","#0090FF","#0094FF","#0098FF","#009CFF","#00A0FF","#00A4FF","#00A8FF","#00ACFF","#00B0FF","#00B4FF","#00B8FF","#00BCFF","#00C0FF","#00C4FF","#00C8FF","#00CCFF","#00D0FF","#00D4FF","#00D8FF","#00DCFF","#00E0FF","#00E4FF","#00E8FF","#00ECFF","#00F0FF","#00F4FF","#00F8FF","#00FCFF","#01FFFD","#05FFF9","#09FFF5","#0DFFF1","#11FFED","#15FFE9","#19FFE5","#1DFFE1","#21FFDD","#25FFD9","#29FFD5","#2DFFD1","#31FFCD","#35FFC9","#39FFC5","#3DFFC1","#41FFBD","#45FFB9","#49FFB5","#4DFFB1","#51FFAD","#55FFA9","#59FFA5","#5DFFA1","#61FF9D","#65FF99","#69FF95","#6DFF91","#71FF8D","#75FF89","#79FF85","#7DFF81","#81FF7D","#85FF79","#89FF75","#8DFF71","#91FF6D","#95FF69","#99FF65","#9DFF61","#A1FF5D","#A5FF59","#A9FF55","#ADFF51","#B1FF4D","#B5FF49","#B9FF45","#BDFF41","#C1FF3D","#C5FF39","#C9FF35","#CDFF31","#D1FF2D","#D5FF29","#D9FF25","#DDFF21","#E1FF1D","#E5FF19","#E9FF15","#EDFF11","#F1FF0D","#F5FF09","#F9FF05","#FDFF01","#FFFC00","#FFF800","#FFF400","#FFF000","#FFEC00","#FFE800","#FFE400","#FFE000","#FFDC00","#FFD800","#FFD400","#FFD000","#FFCC00","#FFC800","#FFC400","#FFC000","#FFBC00","#FFB800","#FFB400","#FFB000","#FFAC00","#FFA800","#FFA400","#FFA000","#FF9C00","#FF9800","#FF9400","#FF9000","#FF8C00","#FF8800","#FF8400","#FF8000","#FF7C00","#FF7800","#FF7400","#FF7000","#FF6C00","#FF6800","#FF6400","#FF6000","#FF5C00","#FF5800","#FF5400","#FF5000","#FF4C00","#FF4800","#FF4400","#FF4000","#FF3C00","#FF3800","#FF3400","#FF3000","#FF2C00","#FF2800","#FF2400","#FF2000","#FF1C00","#FF1800","#FF1400","#FF1000","#FF0C00","#FF0800","#FF0400","#FF0000","#FB0000","#F70000","#F30000","#EF0000","#EB0000","#E70000","#E30000","#DF0000","#DB0000","#D70000","#D30000","#CF0000","#CB0000","#C70000","#C30000","#BF0000","#BB0000","#B70000","#B30000","#AF0000","#AB0000","#A70000","#A30000","#9F0000","#9B0000","#970000","#930000","#8F0000","#8B0000","#870000","#830000","#7F0000"))
options(warn = oldw) 
# silence off

# NUMBERS OF COLORS
clr_num  <-     c(256, 128, 256, 60, 256, 128, 128, 128, 13, 128, 256)


# ABSOLUTE VARIABLES name list
var_names <-     c("interception",  "snowpack",      "SWC_L",         "SM_L",          "SM_Lall",       
                   "sealedSTW",     "unsatSTW",      "satSTW",        "PET",           "aET",           
                   "Q",             "QD",            "QIf",           "QIs",           "QB",            
                   "recharge",      "soil_infil_L",  "aET_L",         "preEffect",     "Qrouted")

# COLOR corresponding to var_names
#                  "cmocean_tempo", "MPL_PuBu","NEO_div_vegetation_a","GMT_drywet",    "GMT_drywet",       
#                  "NCV_blu_red",   "NCV_blu_red",   "NCV_blu_red",   "MPL_RdBu",      "MPL_RdBu",           
#                  "MPL_RdYlBu",    "MPL_RdYlBu",    "MPL_RdYlBu",    "MPL_RdYlBu",    "MPL_RdYlBu",            
#                  "MPL_BrBG",      "precip_diff_12lev","MPL_RdBu",   "MPL_YlGnBu",    "NCV_jet")

# COLOR matrix rows corresponding to var_names
# - lookup vector to link var_names ~ clr_mat, clr_num
clr_mat_row <-   c( 1, 2, 3, 4, 4, 
                    5, 5, 5, 6, 6, 
                    7, 7, 7, 7, 7, 
                    8, 9, 6, 10, 11)

# COLOR DISTRIBUTION corresponding to var_names
# 1 - linear, 3 - more colors for lower values, 0.3 - more colors for higher values
clr_dist      <- c( 3, 0.3, 1, 1, 1,
                    1, 1, 1, 1, 1,
                    0.3, 0.3, 0.3, 0.3, 0.3,
                    3, 3, 1, 3, 3)

# COLOR FLIP corresponding to var_names
clr_flip <-      c("no",            "yes",           "no",            "no",            "no",
                   "no",            "no",            "no",            "yes",           "yes",
                   "yes",           "yes",           "yes",           "yes",           "yes",
                   "no",            "no",            "yes",           "no",            "no")


# SET COLOR for each variable
clr_cutpts = c()
clr_definition = c()

for (ifile in 1:nfiles) { # ================== file loop 
  
  if (nfiles == 2) { # both files available
    if (ifile == 1) {
      ncin <- ncin_mrm
    } else {
      ncin <- ncin_mhm
    }
  } else { # 1/2 files available
    if (exists("ncin_mrm")) {
      ncin <- ncin_mrm
    } else {
      ncin <- ncin_mhm
    }
  }
  
  # list all variables in current file
  varlist <- names(ncin$var)
  varlist <- varlist[ varlist != "time_bnds" & varlist !="lat" & varlist != "lon"] 
  # excluding time and lat lon variables
  
  for (ivar in 1:length(varlist)) { # ========== variable loop
    
    tmp.array <- ncvar_get(ncin,varlist[ivar]) # dimensions (row=lon,col=lat,time) 
    tmp.vector <- as.vector(tmp.array[,,]) # All values (global) of current variable as a vector
    tmp.vector <- na.omit(tmp.vector)
    
    # get suitable COLOR MAP ATTRIBUTES
    idx <- which(var_names == gsub("[!^0-9\\.]", "", varlist[ivar]))
    clr_style <- clr_dist[idx]
    clr_dir <- clr_flip[idx]
    

    # limit the NO. OF CUTPOINTS for the data set to 10000
    if (length(unique(tmp.vector)) <= 10000 ) { # small data set
      ncutpts = length(unique(tmp.vector))-1
    } else {
      ncutpts = 10000
    }
    
    # Find CUT POINTS in the vector of global values &
    # assign and define the color palette
    if (length(unique(tmp.vector)) == 1) { # no simulated value
      intrvls <-  seq(0, 1, length.out =  101)
      length(intrvls) <- 10000
      clr_cutpts <- rbind(clr_cutpts, intrvls)
      clr_map <- colorRampPalette(clr_mat[clr_mat_row[idx],1:clr_num[clr_mat_row[idx]]], bias = clr_style)(100)
    } else {
      # silence on
      oldw <- getOption("warn")
      options(warn = -1) 
        intrvl_object <- classIntervals(tmp.vector, n=ncutpts, style = "equal")
      options(warn = oldw) 
      # silence off
      intrvls <- unique(intrvl_object$brks)
      length(intrvls) <- 10000
      clr_cutpts <- rbind(clr_cutpts, intrvls)
      clr_map <- colorRampPalette(clr_mat[clr_mat_row[idx],1:clr_num[clr_mat_row[idx]]], bias = clr_style)(length(unique(intrvl_object$brks))-1)
      # the bias argument helps to control the COLOR MAP DISTRIBUTION
    }
    
    # FLIP the color?
    if (clr_dir == "yes") {
      clr_map <- rev(clr_map) # reversed
    }
    
    length(clr_map) <- 9999
    clr_definition <- rbind(clr_definition, clr_map)
    # since rbind is row binding in order of variables in the two netCDF  
    # files, while fetching the information, the same order needs to be
    # followed
                        
  } # ========== variable loop close
} # ========== file loop close



#======================================================================== (2) ANIMATE FUNCTION ====

animate <- function(rows, cols, nvar) {
  
  # Progress Bar
  pb=txtProgressBar(min=1, max=nt*nvar, style = 3)
  
  # Initialize
  p <- list() # Defining p as list

  for (itime in 1:nt) { # ======================== time loop
    
    var_count = 0
    
    for (ifile in 1:nfiles) { # ================== file loop 
      
      if (nfiles == 2) { # both files available
        if (ifile == 1) {
          ncin <- ncin_mrm
        } else {
          ncin <- ncin_mhm
        }
      } else { # 1/2 files available
        if (exists("ncin_mrm")) {
          ncin <- ncin_mrm
        } else {
          ncin <- ncin_mhm
        }
      }
      
      # get TIME variable and attributes 
      nctime <- ncvar_get(ncin,"time")
      tunits <- ncatt_get(ncin,"time","units") 
      nt <- dim(nctime)  
      
      tustr <- strsplit(tunits$value, " +")
      tdstr <- strsplit(unlist(tustr)[3], "-")
      tmonth <- as.integer(unlist(tdstr)[2])
      tday <- as.integer(unlist(tdstr)[3])
      tyear <- as.integer(unlist(tdstr)[1])
      tfinal <- chron(nctime/24, origin=c(tmonth, tday, tyear)) 
                # nctime (hours) is converted to days
      
      # Temporal data resolution for date ribbon
      temp_res_check <- tfinal[2] - tfinal[1]
      if (temp_res_check <= 1) {
        date_ribbon <- paste(days(tfinal[itime]), months(tfinal[itime]), years(tfinal[itime]), 
        sep=" ") # date ribbon for daily data
      } else if (temp_res_check >= 32) {
        date_ribbon <- paste(years(tfinal[itime]),sep = " ") 
        # date ribbon for yearly data
      } else {
        date_ribbon <- paste(months(tfinal[itime]), years(tfinal[itime]), sep=" ") 
        # date ribbon for monthly data
      }
      
      # get LAT LON
      
      lon <- ncvar_get(ncin,"lon")
      if (length(dim(lon)) == 2){ # check for 2D lon
        lon <- lon[,1] # we just need a vector of values in the direction in which lon changes
      }
      nlon <- length(lon)
      
      lat <- ncvar_get(ncin,"lat") 
      if (length(dim(lat)) == 2){ # check for 2D lat
        lat <- lat[1,] # we just need a vector of values in the direction in which lat changes
      }
      nlat <- length(lat) 
      
      lat <- rev(lat) # latitude flipping needed for mHM output netCDFs!
      
      # check for Western hemisphere
      if (lon[1] > tail(lon, n=1)) {
        lon <- rev(lon)
      }
      
      # list all variables in current file
      varlist <- names(ncin$var)
      varlist <- varlist[ varlist != "time_bnds" & varlist !="lat" & varlist != "lon"] 
                # excluding time and lat lon variables
      
      
      for (ivar in 1:length(varlist)) { # ========== variable loop
        
        var_count = var_count + 1
        
        # get VARIABLE and its attributes 
        tmp.array <- ncvar_get(ncin,varlist[ivar]) # dimensions (row=lon,col=lat,time) 
        dunits<- ncatt_get(ncin,varlist[ivar],"unit") 
        if (!dunits$hasatt){ # check whether the attribute name for unit is "unit" or "units"
          dunits<- ncatt_get(ncin,varlist[ivar],"units") 
        }
        
        # Prepare grid, color and date
        grid <- expand.grid(lon=lon, lat=lat)
        val <- min(max(lon) - min(lon), max(lon) - min(lon))/2 # Value of suitable interval 
                                                               # (at least 3 cutpoints) along 
                                                               # shorter domain direction
        latlon_interval <- latlon_cutpts[which.min(abs(latlon_cutpts - val))]
        x.scale <- list(at=seq(floor(min(lon)),ceiling(max(lon)),latlon_interval))
        y.scale <- list(at=seq(floor(min(lat)),ceiling(max(lat)),latlon_interval))
        
        
        # PLOT
        tmp.slice <- tmp.array[, , itime]
        tmp.slice <- tmp.slice[,ncol(tmp.slice):1] # flip the latitude values as "lat" was flipped
        myplot <- levelplot(tmp.slice ~ lon * lat, data=grid, aspect = "iso", pretty=T, 
                            at=na.omit(clr_cutpts[var_count,]), col.regions=na.omit(clr_definition[var_count,]), 
                            xlab = "", ylab = "",
                            main = list(paste(varlist[ivar]," (",dunits$value,")",sep=""), 
                                        cex=1.4 - 0.05*(rows-1)), 
                            colorkey = list(tick.number=5, labels=list(cex=1.4 - 0.1*(rows-1))), 
                            scales=list(tck=c(-1,-1), x=list(x.scale, cex=1.2 - 0.1*(rows-1)), 
                                        y=list(y.scale, cex=1.2 - 0.1*(rows-1)), xlab="", ylab=""))
        
        # Save plot as an element of multiplot
        xplot = gridExtra:::latticeGrob(myplot)
        p[[ivar + (ifile-1) ]] <- xplot # p[[1]] is dedicated to the empty plot
        
        # Progress bar update
        setTxtProgressBar(pb, (itime-1)*nvar + ivar)
        
        
      } # ========== variable loop close
    } # ========== file loop close

    
    # Determine the number of training blank plots for the multiplot
    count_empty_plots <- cols * rows - nvar    # number of empty plot to pad to the layout
    
    if (count_empty_plots != 0){
      xlay <- matrix(c(seq(1,length(p[])),replicate(count_empty_plots,NA)), nrow = rows, ncol = cols, 
        byrow = TRUE)
    } else { # 'replicate' creates problem when count_empty_plots = 0
      xlay <- matrix(c(seq(1,length(p[]))), nrow = rows, ncol = cols, byrow = TRUE)
    }
    
    select_grobs <- function(lay) {
      id <- unique(c(lay)) 
      id[!is.na(id)]
    }
    
    mymultiplot <- grid.arrange(grobs=p[select_grobs(xlay)], layout_matrix=xlay, 
                                top = textGrob(paste("\n", date_ribbon,"\n", sep=" "),
                                               gp=gpar(fontsize=22 - 0.5*(rows-1),font="Helvetica", 
                                                just="left")))
  
  } # ========== time loop close
  close(pb)
} # ============ function end



#============================================================================= (3) OUTPUT =========

# (3.1) PDF generation ==========================

oldw <- getOption("warn")
options(warn = -1)

pdf( paste("animate.pdf", sep=""), width = 9 + 2*(cols-1), height =  6 + 2*(rows-1))
print("generating pdf...")

animate(rows, cols, nvar) # call animate

options(warn = oldw)



# (3.2) GIF generation ==========================

## ! NOTE: make sure ImageMagick has been installed in your system

oldw <- getOption("warn")
options(warn = -1)

saveGIF({  
  print("generating gif...")

  animate(rows, cols, nvar) # call animate

}, movie.name = paste("animate.gif", sep=""), interval=gifinterval, nmax = 500, 
ani.width = 600 + 200*(cols-1), ani.height = 400 + 200*(rows-1))

options(warn = oldw)


