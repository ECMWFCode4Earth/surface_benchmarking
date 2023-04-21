Introduction:

Welcome to the surface benchmarking challenge, which uses the ECMWF land surface validation tool called "landver" (Fairbairn et al., 2019). An example is setup here to validate soil moisture and soil temperature over 2020 using in situ observations from the international soil moisture network (Dorigo et al., 2011). In situ data from the SMOSMANIA network has also been provided directly by Meteo-France (Calvet et al., 2008). Two analysis products are considered here (i) the ERA5 reanalysis (labelled 0001_ea) (Hersbach et al., 2020), which is freely available on the Copernicus climate data store (C3S, 2023); (ii) research experiment with an updated vegetation mapping (labelled hyfs_rd). For each of these experiments, soil moisture and soil temperature files are provided at 00 and 12 UTC over 2020. Below you will find instructions on how to carry out the SM/ST validation for these datasets. Additional Information and references for the in situ network providers can be found in the Section "In situ data providers" at the end of this document. Additional instructions on landver and interpreting the results can be found here: https://confluence.ecmwf.int/display/~dadf/LANDVER

Instructions:

1. Setup a local conda environment and install the necessary software by downloading and running one of the following scripts in this github repository (i) setup_miniconda_linux (for linux users) or (ii) setup_miniconda_mac (for mac users). Unfortunately the installation for windows is not currently available. Note that each time you run landver you need to activate the local environment by following steps 3 and 5.    

2. Download the data (both sfc_evaluation_out.zip, containing the analysis SM/ST data and ISMN_data.zip, containing the in situ SM/ST data) from the following link (https://myftp.ecmwf.int/files/public/land_surface_data/) and unzip both folders. 

3. Download the landver python code on this github repository (landver_code.zip) and unzip. Following the instructions in the README file in the landver_code folder to run the landver example validating SM and ST. The results can be found in the html output. Some examples intepreting other experiments can be found on the confluence page (https://confluence.ecmwf.int/display/~dadf/LANDVER) or in Fairbairn et al. (2019). 


In situ data providers:

The analysis soil moisture is validated against in situ observations form the ISMN (Dorigo et al., 2011), which are already extracted and preprocessed in the insitu_dir in LDAS_config. The in situ observations and potentially erroneous observations are screened (due to frozen conditions, spurious spikes, etc...). The networks consist of SCAN, USCRN and SNOTEL in the US, SMOSMANIA in France, REMEDHUS in Spain and Oznet in Australia. 

    SCAN (Schaefer et al., 2007), US (up to 184 stations): Obs depths (10, 20, 50, 100 cm)
    
    USCRN (Bell et al., 2013), US (up to 109 stations): Obs depths (10, 20, 50, 100 cm)
    
    SNOTEL (Schaefer et al., 2007), US (up to 380 stations): Obs depths (10, 20, 50, 100 cm)
    
    SMOSMANIA (Calvet et al., 2007, Albergel et al., 2008), France (up to 21 stations): Obs depths (0, 10, 20, 30 cm)
    
    REMEDHUS (Martínez-Fernández and Ceballos, 2005), Spain (Up to 16 stations): Obs depths (5 cm)
    
    OZNET (Smith et al., 2012), Australia (Up to 37 stations): Obs depths (4, 15,45, 75 cm)
    
    TERENO (Zacharias et al., 2011, Bogena et al., 2018), Germany (Up to 4 stations): Obs depths (5, 20, 50 cm)

References:

C. Albergel, C. Rudiger, T. Pellarin, J.-C. Calvet, N. Fritz, F. Froissard, D. Suquia, A. Petitpa, B. Piguet, and E. Martin, “From near-surface to root-zone soil moisture using anexponential filter: an assessment of the method based on in situ observations and model simulations,” Hydrol. Earth Syst. Sci., vol. 12, pp. 1323–1337, 2008.

J. E. Bell, M. A. Palecki, C. Baker, W. Collins, J. Lawrimore, R. Leeper, M. Hall, J. Kochendorfer, T. Meyers, T. Wilson, and H. Diamond, “U.S. Climate Reference Network soil moisture and temperature observations,” J. Hydrometeor., vol. 14, pp. 977–988, 2013.

Bogena, H.R., Montzka, C., Huisman, J.A., Graf, A., Schmidt, M., Stockinger, M., Von Hebel, C., Hendricks-Franssen, H.J., Van der Kruk, J., Tappe, W. and Lücke, A., 2018. The TERENO‐Rur hydrological observatory: A multiscale multi‐compartment research platform for the advancement of hydrological science. Vadose Zone Journal, 17(1), pp.1-22, https://doi.org/10.2136/vzj2018.03.0055;

J. Calvet, N. Fritz, F. Froissard, D. Suquia, A. Petitpa, and B. Piguet, “In situ soil moisture observations for the CAL/VAL of SMOS: the SMOSMANIA network.” in Geoscience and Remote Sensing Symposium, IGARSS 2007., vol. 16 (3). IEEE International, 2007, pp. 1293–1314.

Copernicus Climate Change Service (C3S) (2017): ERA5: Fifth generation of ECMWF atmospheric reanalyses of the global climate . Copernicus Climate Change Service Climate Data Store (CDS), Last accessed April 2023. https://cds.climate.copernicus.eu/cdsapp#!/home

Dorigo, W.A., Wagner, W., Hohensinn, R., Hahn, S., Paulik, C., Xaver, A., Gruber, A., Drusch, M., Mecklenburg, S., van Oevelen, P. and Robock, A., 2011. The International Soil Moisture Network: a data hosting facility for global in situ soil moisture measurements. Hydrology and Earth System Sciences, 15(5), pp.1675-1698.

Fairbairn, D., de Rosnay, P., & Browne, P. A. (2019). The new stand-alone surface analysis at ECMWF: Implications for land–atmosphere DA coupling. Journal of Hydrometeorology, 20(10), 2023-2042. doi: https://doi.org/10.1175/JHM-D-19-0074.1 

Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz‐Sabater, J., Nicolas, J., Peubey, C., Radu, R., Schepers, D. and Simmons, A., 2020. The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), pp.1999-2049. doi: https://doi.org/10.1002/qj.3803

J. Martínez-Fernández and A. Ceballos, “Mean soil moisture estimation using temporal stability analysis,” J. Hydrol., vol. 312, pp. 28–38, 2005.

G. Schaefer, M. Cosh, and T. Jackson, “The USDA natural resources conservation service soil climate analysis network (SCAN),” J. Atmos. Oceanic Technol., vol. 24(2), pp. 2073–2077, 2007.

Smith, A. B., Walker, J. P., Western, A. W., Young, R. I., Ellett, K. M., Pipunic, R. C., Grayson, R. B., Siriwidena, L., Chiew, F. H. S. and Richter, H. The Murrumbidgee Soil Moisture Monitoring Network Data Set. Water Resources Research, vol. 48, W07701, 6pp., 2012 doi:10.1029/2012WR011976 

Zacharias, S., H.R. Bogena, L. Samaniego, M. Mauder, R. Fuß, T. Pütz, M. Frenzel, M. Schwank, C. Baessler, K. Butterbach-Bahl, O. Bens, E. Borg, A. Brauer, P. Dietrich, I. Hajnsek, G. Helle, R. Kiese, H. Kunstmann, S. Klotz, J.C. Munch, H. Papen, E. Priesack, H. P. Schmid, R. Steinbrecher, U. Rosenbaum, G. Teutsch, H. Vereecken. 2011. A Network of Terrestrial Environmental Observatories in Germany. Vadose Zone J. 10. 955–973, https://doi.org/10.2136/vzj2010.0139;

