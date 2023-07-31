# Benchmarking Surface Heat Fluxes

Validation of soil moisture and soil temperature is crucial for Numerical Weather Prediction (NWP), as they control surface heat fluxes that directly affect near-surface weather. This can be done with LANDVER, which is a validation package for land surface variables, currently consisting of soil moisture and soil temperature. The tool provides an independent validation of soil moisture and soil temperature data using in situ observations from the international soil moisture network (Dorigo et al., 2011). What is currently missing in this software package is the capability to validate latent and sensible surface heat fluxes against Eddy-Covariance measurements, which can provide useful information about how well ECMWF’s Land-Surface Modelling Component ECLand is able to translate soil moisture stress into surface heat fluxes. Implementing an additional routine into the already existing software package paves the way for a standardized land-surface benchmarking tool for the ECMWF.

Details: https://github.com/ECMWFCode4Earth/challenges_2023/issues/12

