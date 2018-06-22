# Model interface GIS2FEM #

Copyright: Wenqing Wang (wenqing.wang@ufz.de)
           Miao Jing (miao.jing@ufz.de)

Online Repository: https://github.com/UFZ-MJ/OGS_mHM

Citation: Jing, M., Heße, F., Kumar, R., Wang, W., Fischer, T.,
Walther, M., Zink, M., Zech, A., Samaniego, L., Kolditz, O.,
and Attinger, S.: Improved regional-scale groundwater
representation by the coupling of the mesoscale Hydrologic
Model (mHM v5.7) to the groundwater model OpenGeoSys (OGS),
Geosci. Model Dev., 11, 1989-2007,
https://doi.org/10.5194/gmd-11-1989-2018, 2018. 

GIS2FEM has been embedded into the source code of OGS5 (https://www.opengeosys.org/).
To get the excutable file of GIS2FEM, please download also the source code of OGS5 via the following link:
https://github.com/UFZ-MJ/OGS_mHM/tree/master/OGS5.
After compling the source code,
the executable file of GIS2FEM, as well as the executable file of OGS5 can be obtained.
 
Tutorial:
This is a brief instruction of how to convert the triangular-wise or quadrilateral-wise recharge data 
into the nodal source terms of a three dimensional finite element model for OGS.

Assume the file name of the  the 3D mesh file name is 'Naegelstedt.msh', 
and the file containing the mHM-generated raster file of the gridded recharge is 
'recharge-monthly.asc'.

The data process takes the following procedure.
      
First, a new input file with the suffix of "pcp" has to be prepared. The file has four lines as:

	Unit: Month
	Ratio: 1
	recharge-monthly.asc

The first line have a keyword 'Unit' and it is for time unit.
The second line has one keyword 'Ratio' and one value. The value can be used to 
scale the recharge data. For instance, the unit of the recharge data is mm/d 
but the unit of m/d is required in the simulation by OGS, then we set 1.e-3 
for 'Ratio'. The third line gives the name of the mHM-generated raster file of the gridded recharge.
Assuming this new file is saved as 'monthly.txt', the data 
conversion can then be performed as the last step.

Second, with the .pcp file prepared in the previous step, the recharge data 
is then converted into the OGS input files by using the command with option '1' as 
   GIS2FEM 1 Naeglestedt
Program GIS2FEM first reads the top surface elements, which could be triangle or 
quadrilateral, and the gridded recharge data. Then it assign the proper 
recharge value to the top surface elements by checking the coordinates of the centroid 
of each top surface element. For each surface element, if its centroid is in a triangle 
of the grid of  the recharge data, the value of the triangle is assigned to the surface 
element.  Once all top surface elements have been processed, the elements that have been 
assigned with the recharge values are involved in the face integration calculation, 
in which the recharge data is converted into nodal source terms. The calculated results 
are written to a binary file called "recharge-monthly.asc.bin".
This file can be directly read by OGS and served as Neumann BC of the OGS model.
 Besides this file,  one more control file  is generated as well after the command is done. 
The additional file is named automatically by adding string '.ifl' to the name of 
the project. In this example, it is
   Naegelstedt.ifl
The file starts with the same header of the raster file of the gridded
recharge, and follows with the names of the generated 
files and their corresponding times, e.g.
	
	ncols 80
	nrows 72
	xllcorner 4.374e+06
	yllcorner 5.656e+06
	cellsize 500
	NODATA_value -9999
	0 recharge-monthly.asc.bin
	#STOP

OGS simulation:
For the simulation, one only needs to fill '$DIS_TYPE' in .st input
file with key word 'PRECIPITATION'  and the name of the control file
name that generated  the last step of the above data conversion. 

Here is an example:
 
	#SOURCE_TERM
	$PCS_TYPE
	GROUNDWATER_FLOW
	$PRIMARY_VARIABLE
	HEAD
	$DIS_TYPE
	PRECIPITATION Naegelstedt.ifl
	#STOP
