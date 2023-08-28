# Bayesian statistical modeling to predict observer-specific optimal windowing parameters in magnetic resonance imaging

This repository contains the source code and parameters for the Bayesian statistical model used in the research paper titled **"Bayesian statistical modeling to predict observer-specific optimal windowing parameters in magnetic resonance imaging**" by Kohei Sugimoto, Masataka Oita, and Masahiro Kuroda, published in *Heliyon*  2023 (doi: [10.1016/j.heliyon.2023.e19038](https://doi.org/10.1016/j.heliyon.2023.e19038)).


## Method overview

### Training

1. Prepare the dataset of raw histograms from MRI images for each protocol and MRI operator.
2. Create mapping functions from raw histograms to standard histograms for each protocol, applying Nyul normalization.
3. Map windowing parameters to standard histograms for each protocol to obtain standardized windowing parameters.
4. Fit the standardized windowing parameters to a Bayesian statistical model and obtain the optimal standardized windowing parameters for each protocol and MRI operator.

### Prediction

1. Prepare the MRI image for which you want to obtain the optimal windowing parameters.
2. Specify the protocol and select an arbitrary MRI operator (who has set the windowing parameters).
3. Map the standardized optimal windowing parameters from the protocol and MRI operator selected in step 2 to the MRI images prepared in step 1.
4. Obtain the optimal windowing parameters at the original scale using the inverse function of the mapping function.

## Data Locations

- The Bayesian statistical model use in the paper, written in *Stan*, is located in `model/model.stan`.
- Standardized histograms generated using Nyul's standardizization can be found in `parameters/standardized_histograms`.
- Standardized windowing parameters are located in `parameters/standardized_windowing_parameters.csv`.

## Requirements

- Programming language
    - Python 3.7.3
- Packages
    - numpy
    - pydicom
    - nibabel 
    - pystan 
    - intensity-normalization

