# LANDVER for surface fluxes


### Scripts
├── `Common_functions.py`: class LDAS_config, that contains default parameters; functions for html output; parser <br>
├── `landver`: main <br>
├── `LANDVER_pre.py`: preprocessing of data<br>
├── `LANDVER.py`: in-situ validation<br>
├── `Options_namelist.py`: user-defined variables, e.g. experiments, analysis period<br>
└── `Validation_functions.py`: utilities for validation and plotting


### How does the code work?
- namelist parameters are read in (namelist: `Options_namelist.py`) and converted (with `getConfigutation` from `Common_functions.py`) into an instance of the class `LDAS` (land data assimilation) called `cfg`
- preprocessing of model analyses (here: 0001_ea, hyfs_rd) in function `preprocessData` in `LANDVER_pre.py`: reads model data from grib files, merges them for the desired analysis period, for each station (from station metadata file) the nearest model grid point is found, conversion of units, deaccumulation -> per station a new output file (.dat) is created, which contains date, sensible and latent heat flux (after all preprocessing in W/m^2)
- Validation of modelled fluxes against measured fluxes in function `in_situ_validation` in `LANDVER.py`:
    1. initialize dictonaries to store scores
    2. calculate anomaly correlation coefficient (ACC) based on 15-day moving window
    3. calculate scores: loop over experiments, loop over network -> scores: bias, rmse/rmsd, unbiased rmsd, correlation, acc
        - over entire period
        - check significance (if `cfg.stat_sig`) based on p-value
        - for every season (from `slist`) -> seasonal scores
        - for every year (in `years.year`) -> annual scores
        - plot time series (if `cfg.plot_time_series`) with `plot_time_series` function from `Validation_functions.py`
        - filter scores based on p-value
        - update dictionaries
    4. check that all experiments are validated for the same stations
    5. generate tables and boxplots (of the previously calculated scores) for each network and experiment separately
        - station average
        - check whether experiment performs (significantly) better than control experiment
        - create directories to store results
    6. calculate scores for all networks together (average over all networks)
- create html file that summarizes all scores (with function `web_create` from `Common_functions.py`)



├── `Common_functions.py` <br>
├── `landver`<br>
├── `LANDVER_pre.py`<br>
├── `LANDVER.py`<br>
├── `Options_namelist.py` <br>
└── `Validation_functions.py`
