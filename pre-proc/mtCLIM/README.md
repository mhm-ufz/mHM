## mtCLIM preprocessor 
## version 1.001

Adapted pre-processing of meteorological input for the mesoscale hydrological model. 
This Fortran module is based on [mtCLIM v4.3](http://www.ntsg.umt.edu/project/mt-clim.php)  

Taking into account base meteorological variables (minimum and maximum air temperature, precipitation) and 
morphological characteristics of the underlying terrain (digital elevation model, slope, aspect) 
the program is able to estimate humidity (vapore pressure or vapore pressure deficit) and incoming shortwave radiation. 

Allowing for albedo as an additional input variable it extends mtCLIM v4.3 by the estimation of net shortwave and longwave radiation. 
Thus, providing an estimate of net radiation (and actual vapore pressure) as important mHM input variable, 
e.g. for radiation based approches to estimate potential evapotranspiration.  

Albedo data has to be given as a 12-monthly climatology at L0 level (see gridded LAI data) including white and black sky albedo
(source: e.g [MODIS](https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mcd43a3))

#### how to start?
Have a look to the _ini_ namelist. Adjust file paths according to your needs.  
Compile *f90 code by simply running `make`. Run mtclim `./mtclim_bin`.

#### dependencies 
FORTRAN modules in FORTRAN_chs_lib
* mo_append.f90
* mo_constants.f90
* mo_julian.f90
* mo_kind.f90
* mo_message.f90
* mo_mtclim.f90
* mo_ncread.f90
* mo_netcdf.f90
* mo_string_utils.f90
* mo_utils.f90  

#### contact
If any problems arise with this piece of code, please send me a message!  
Thanks.  
Johannes Brenner, CHS - UFZ Leipzig, johannes.brenner@ufz.de  

current version according to last gitlab commit:  
c9aaae4c, origin master, 2018-01-19
