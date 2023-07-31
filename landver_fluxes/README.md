# LANDVER for surface fluxes


### scripts
- run: `./landver Options_namelist.py`
- `Options_namelist.py`: user-defined variables, e.g. analysis period
- `Common_functions.py`: class LDAS_config -> default parameters; contains functions for html output
- `landver`: main
- `LANDVER_pre.py`: preprocessing of data
- `LANDVER.py`: in-situ validation function
- `Validation_functions.py`: contains the actual validation functions; functions for reading the data; plot functions


### How does the code work?
- namelist parameters are read in (namelist: `Options_namelist.py`) and converted (with `getConfigutation` from `Common_functions.py`) into an instance of the class `LDAS` (land data assimilation) called `cfg`
- preprocessing of model analyses (here: ea_0001, hyfs_rd) in function `preprocessData` in `LANDVER_pre.py`: reads model data from grib files, merges them for the desired analysis period, for each station (from station metadata file) the nearest model grid point is found, conversion of units, deaccumulation -> per station a new output file (.dat) is created, which contains date, sensible and latent heat flux (after all preprocessing in W/m^2)
- Validation of modelled fluxes against measured fluxes in function `in_situ_validation` in `LANDVER.py`