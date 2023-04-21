Welcome to the surface benchmarking challenge, which uses the ECMWF land surface validation tool called "landver" (Fairbairn et al., 2019). An example is setup here to validate soil moisture and soil temperature over 2020 using in situ observations from the international soil moisture network (Dorigo et al., 2011). In situ data from the SMOSMANIA network has also been provided directly by Meteo-France (Calvet et al., 2008). Two analysis products are considered here (i) the ERA5 reanalysis (labelled 0001_ea) (Hersbach et al., 2020), which is freely available on the Copernicus climate data store (C3S, 2023); (ii) research experiment with an updated vegetation mapping is provided (labelled hyfs_rd). For each of these experiments, soil moisture and soil temperature files are provided at 00 and 12 UTC over 2020. Below you will find instructions on how to carry out the SM/ST validation for these datasets. Additional Information and references for the in situ network providers can be found in the Section "In situ data providers" at the end of this document. Additional instructions on landver and interpreting the results can be found here: https://confluence.ecmwf.int/display/~dadf/LANDVER

Instructions:

1. Setup a local conda environment and install the necessary software by downloading and running one of the following scripts in this github repository (i) setup_miniconda_linux (for linux users) or (ii) setup_miniconda_mac (for mac users). Unfortunately the installation for windows is not currently available. Note that each time you run landver you need to activate the local environment by following steps 3 and 5.    

2. Download the data (both sfc_evaluation_out.zip, containing the analysis SM/ST data and ISMN_data.zip, containing the in situ SM/ST data) from the following link (https://myftp.ecmwf.int/files/public/land_surface_data/) and unzip both folders. 

3. Download the landver python code on this github repository (landver_code.zip) and unzip. Following the instructions in the README file in the landver_code folder to run the landver example validating SM and ST. The results can be found in the html output. Some examples intepreting other experiments can be found on the confluence page (https://confluence.ecmwf.int/display/~dadf/LANDVER) or in Fairbairn et al. (2019). 


In situ data providers:

The analysis soil moisture is validated against in situ observations form the ISMN (Dorigo et al., 2011), which are already extracted and preprocessed in the insitu_dir in LDAS_config. The in situ observations and potentially erroneous observations are screened (due to frozen conditions, spurious spikes, etc...). The networks consist of SCAN, USCRN and SNOTEL in the US, SMOSMANIA in France, REMEDHUS in Spain and Oznet in Australia. 
    SCAN, US (up to 204 stations): Obs depths (0,10,20,50,100 cm)
    USCRN, US (up to 78 stations): Obs depths (0,10,20,50,100 cm)
    SNOTEL, US (up to 418 stations): Obs depths (0,10,20,50,100 cm)
    SMOSMANIA, France (up to 21 stations): Obs depths (0,10,20, 30 cm)
    REMEDHUS, Spain (Up to 24 stations): Obs depths (0.05 cm)
    OZNET, Australia (Up to 37 stations): Obs depths (0.04,0.15,0.45, 0.75 cm)
    TERENO, Germany (Up to 4 stations): 2019-2021

References:

Copernicus Climate Change Service (C3S) (2017): ERA5: Fifth generation of ECMWF atmospheric reanalyses of the global climate . Copernicus Climate Change Service Climate Data Store (CDS), Last accessed April 2023. https://cds.climate.copernicus.eu/cdsapp#!/home

Dorigo, W.A., Wagner, W., Hohensinn, R., Hahn, S., Paulik, C., Xaver, A., Gruber, A., Drusch, M., Mecklenburg, S., van Oevelen, P. and Robock, A., 2011. The International Soil Moisture Network: a data hosting facility for global in situ soil moisture measurements. Hydrology and Earth System Sciences, 15(5), pp.1675-1698.

Fairbairn, D., de Rosnay, P., & Browne, P. A. (2019). The new stand-alone surface analysis at ECMWF: Implications for land–atmosphere DA coupling. Journal of Hydrometeorology, 20(10), 2023-2042. doi: https://doi.org/10.1175/JHM-D-19-0074.1 

Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz‐Sabater, J., Nicolas, J., Peubey, C., Radu, R., Schepers, D. and Simmons, A., 2020. The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), pp.1999-2049. doi: https://doi.org/10.1002/qj.3803
