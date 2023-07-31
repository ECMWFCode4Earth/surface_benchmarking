# LANDVER for surface fluxes


### scripts
- Options_namelist.py: user-defined variables, e.g. analysis period
- Common_functions.py: class LDAS_config -> default parameters; contains functions for html output
- run: ./landver Options_namelist.py
- landver: main
- LANDVER_pre.py: preprocessing of data
- LANDVER.py: in-situ validation function
- Validation_functions.py: contains the actual validation functions; functions for reading the data; plot functions
