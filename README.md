#lake-ice-albedo
Google Earth Engine and Matlab scripts used for tracking lake and land snow fraction and albedo and lake contribution to the terrestrial Crysophere Radiative Effect. Any questions about these scripts should be directed to project leader Sarah Cooley at sarah.cooley@duke.edu 

### **Step 1: Process MODIS data in Google Earth Engine**

Script: run\_lake\_tiles\_dec17\_v1.ipynb (implemented in Google CoLab)

Description: This script runs all the processing of MODIS data (both snow cover and albedo) for each 1° x 1° grid cell. Requires inputs of grid cell shapefiles and separate eroded lake (for lake pixels) and buffered lake (to exclude for land pixels) shapefiles. This script outputs separate lake and land time series csvs for each grid cell for each year, each containing 7 variables, namely daily (1) water pixel count, (2) snow pixel count, (3) cloud pixel count, (4) snow albedo, (5) water albedo, (6) snow albedo standard deviation, (7) water albedo standard deviation.  

### **Step 2: Read in CSV files and save into Matlab format**

Run Script: run\_read\_in\_csv\_time\_series\_v1.m

Contained Functions: read\_in\_csv\_time\_series\_v1.m

Description: his script reads in the GEE-exported CSV files, fill in any missing days with NaN values, and check for any missing years/days in the resulting time series. 

### **Step 3: Process time series**

Run script: process\_results\_monthly\_sep2.m

Contained Functions:

* filter\_ts.m
* filter\_ts_albedo.m
* get\_monthly\_values.m

Description: This script converts the raw daily time series into processed monthly data. This includes removing outliers and cloud-obscured observations, filtering and interpolating to daily, calculating monthly values, and converting to radiative effect through combining results with kernel data. 

### **Step 4: Complete analysis and make figures**

Run script: complete\_processing\_and\_analysis\_sep2.m

Contained Functions:

* calculate\_statistics.m 
* calculate\_statistics\_monthly.m
* get\_global\_weighted\_lc\_stats.m
* get\_global\_weighted\_stats.m
* get\_global\_weighted\_lakecont\_stats.m
* get\_global\_weighted\_lakeinc\_stats.m
* get\_global\_stats\_tables.m
* get\_biome\_stats\_tables.m
* get\_sensitivity\_table.m
* calculate\_sensitivity.m
* variables\_for\_fillplot\_v2.m
* make\_map\_plot.m

Description: This script completes the analysis for all results presented in the manuscript. This includes calculating all statistics for each grid cell, for the whole northern hemisphere domain and by biome, converting from radiative effect to radiative forcing by calculating the slope of the relationship between air temperature and radiative effect, and creating all figures.  
