# Benchmarking Surface Heat Fluxes

**Project Description**: Validation of soil moisture and soil temperature is crucial for Numerical Weather Prediction (NWP), as they control surface heat fluxes that directly affect near-surface weather. This can be done with LANDVER, which is a validation package for land surface variables, currently consisting of soil moisture and soil temperature. The tool provides an independent validation of soil moisture and soil temperature data using in situ observations from the international soil moisture network (Dorigo et al., 2011). What is currently missing in this software package is the capability to validate latent and sensible surface heat fluxes against Eddy-Covariance measurements, which can provide useful information about how well ECMWF's Land-Surface Modelling Component ECLand is able to translate soil moisture stress into surface heat fluxes. Implementing an additional routine into the already existing software package paves the way for a standardized land-surface benchmarking tool for the ECMWF.

**Aim**: Development of a Python package ("LANDVER") for validating ECMWF's surface modules against insitu measurements

**Details**: https://github.com/ECMWFCode4Earth/challenges_2023/issues/12

## Repository Structure
- landver_code.zip: original version of landver for soil temperature and soil moisture
- landver_soil: contains unzipped python code of landver for validation of soil temperature and soil moisture
- landver_fluxes: adapted landver code for validation of surface heat fluxes (sensible and latent heat flux)
- notebooks: contains a collection of jupyter notebooks for analyzing and validating surface heat fluxes against eddy-covariance measurements from ICOS

## Data

#### Model Data:
- **ERA5** (control): coupled land-atmosphere model with data assimilation, climate: v015, resolution: ~31 km, 1 hour
- **ERA5-Land** (experiment): cycle: 45r1(2018), resolution: ~ 9 km, 1 hour
- **HyFS** (experiment): offline surface model run, cycle: 49r1 (new), climate: v021, resolution: ~ 31 km, 6 hours

#### Insitu-Data:
- **ICOS** (Integrated Carbon Observation System, https://www.icos-cp.eu/): eddy-covariance measurements from flux towers in Europe (point measurements, temporal resolution: 30 min)


## Installation
1. Clone the github repository with `git clone git@github.com:ECMWFCode4Earth/surface_benchmarking.git`
2. Set up the environment by running `setup_miniconda_linux`
3. Activate the environment with
    ```
    export PATH="$HOME/miniconda/bin:$PATH"
    source activate $HOME/miniconda/envs/landver_env/
    ```
4. Run LANDVER or one of the jupyter notebooks (detailed explanation in the respective folders)

