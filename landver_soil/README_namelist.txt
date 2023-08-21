Welcome to the land verification tool, which was demonstrated by Fairbairn et al., (2019). The tool provides an independent validation of soil moisture and soil temperature data using in situ observtions from the international soil moisture network. More detailed information is also available 
on the confluence page: 
https://confluence.ecmwf.int/display/~dadf/LANDVER 

The following instructions are for running landver with the ECMWFcode4Earth challenge 15:


1. Extracting and preparing the soil moisture analysis data:

    1.1 The required user-defined variables are in the file Options_namelist.py.  
        Select the period you wish to extract/validate the SM/ST analysis data:

        analysis_period = ['2020-01-01','2020-12-31] #Analysis period.
 

    1.2 The Options_namelist.py will require available data over the analysis period. Change the variable output_dir to point to the directory where the analysis data is stored         i.e. sfc_evaluation_out  

        The experiments can have different configurations, but must contain daily global soil moisture and soil temperature outputs.
        Make sure all the fields are filled in for each of the experiments in "EXPVER". The following example validates ERA5 and the RD experiment hyfs. 
        The first experiment (in this case ERA5) is the "control".
    

	#Experimental parameters
	EXPVER =    ["0001","hyfs"] #Experiment names needed for extraction, preprocessing and validation 
	#First experiment is considered the control 
	CLASS =     ["ea","rd"]
	TYPE =      ["AN","AN"]
	ANOFFSET =  ["9","9"]
	STREAM =    ["oper","oper"]
	REPRES =    ["SH","SH"]
	LEVTYPE =   ["SFC","SFC"]
	TIME      = ["00/12","00/12"] 
	STEP =      ["00","00"]
	DOMAIN =    ["G","G"]
	GRID =      ["av","av"]
	EXPNAME =  ["ERA5","hyfs"]


    1.3 Check the default parameters in the file Common functions.py, Class LDAS_config. Change the parameters, if preferred. 

    - In Common_functions.py, change the directory insitu_dir to the folder with the ISMN observations, i.e. .../ISMN_data/in_situ_data/pre_processed_data

    - Both soil moisture and soil temperature raw analysis data are already extracted for this example. They will be preprocessed and validated using the following options:

    self.extract_SM = False  # Extract soil moisture grib files from MARS (required for preprocessing)
    self.extract_ST = False 
    self.pre_process_SM = (
        True  # Preprocess soil moisture (required for validation)
    )
    self.pre_process_ST = (
        True  
    )
    self.validate_SM = True  # Validate soil moisture
    self.validate_ST = True  

    - to avoid preprocessing the analysis data twice, the self.pre_process_SM/self.pre_process_ST can be set to False after running landver for the first time.

    1.4 Select the in situ networks where you want to preprocess the analysis data (nearest neighbour) e.g.
    Network = ['USCRN', 'SNOTEL', 'SCAN', 'SMOSMANIA', 'OZNET','REMEDHUS','TORENO']  # Can be any combination of 'USCRN', 'SNOTEL', 'SCAN', 'SMOSMANIA', 'REMEDHUS', 'OZNET',
               'TERENO'

    Please see the confluence page for information on the other parameters.

2.  Run the preprocessing/validation: 
    - Make sure your local environment is active with the python environment:
	
        export PATH="$HOME/miniconda/bin:$PATH"
	source activate landver_env

    - Perform the preprocessing/validation:
    ./landver Options_namelist.py

3.  At the end of the validation an html output is generated, which includes plots and tables of the scores (Pearson R, RMSD, unbiased RMSD, bias)
    See the confluence page above for an example.












