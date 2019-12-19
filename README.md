# Wavelet analysis of directional data

Example of how to perform wavelet analysis on a sequence of directional data, using the example of the flow direction of the Colorado river, USA. Supplementary code relating to the paper *"Flow directions of drainage networks set by the mantle"* (sub judice) 

## Requisites

DEM preprocessing was written for [Generic Mapping Tools](http://gmt.soest.hawaii.edu/) v6 but is compatible with some older versions. 

The python scripts require `python` version **3.6.9** and a modified version of the `mlpy` package. Generic scientific computing pacakges are also used (e.g. `numpy`, `pandas`, `matplotlib` etc...)

## Contents

This repository contains 5 files. 

### Data preprocessing

`preprocessing.sh` indicates how to resample a river stored as a sequence of longitudes and latitudes into a sequence of azimuth as a function of downstream distance. This `bash` script calls a `python` script, 

`find_near.py` is an auxilliary script called from within `preprocessing.sh` which finds the nearest adjacent point for a given river.

### Wavelet Analysis

`azimuths_wavelets_colorado.py` performs wavelet transformation on the Colorado river flow-direction as a function of distance downstream by transforming azimuth into the complex plane. The script transforms the bearing into a function of distance and wavenumber, and subsequently displays the results. The script then performs the inverse wavelet transformation to reconstruct the original signal, as well as reconstructing filtered versions of the signal by removing variation contained on small wavelengths. Finally, the script shows that performing wavelet analysis on azimuths by splitting it into 'eastings' and 'northings' produces identical results to the method in the complex plane.   

`colorado.dat` is a text file containing the longitude, latitude, distance downstream and flow direction of the Colorado river. 

`continuous.py` is a **modified** source code of the `mlpy` package. These modifications can be implemented by simply replacing the original `continuous.py` (found where your python implementation is stored) with the one in this repository. Alternatively, the modifications are minor and are detailed in `continuous.py`.

