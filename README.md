# calculate_Ejk
Overview

This repository contains the core Python code for calculating (E_{jk}), the root-mean-square (RMS) differences between weather forecasts initialized on different days but verifying on the same valid day. The calculation is used to characterize forecast error growth and assess the predictability of operational weather forecasting systems.

The methodology and its application to operational forecasting systems are described in the associated paper.

Data

Users must prepare the input data independently. This repository does not include the forecast datasets used in the study.

To obtain the relevant publicly available data, please refer to the data sources and descriptions provided in the paper.

The input data should contain the geopotential height at 500 hPa (z500) forecasts required for the (E_{jk}) calculation. Please ensure that the forecast initialization dates, lead times, valid dates, spatial grids, and units are consistent with the requirements of the code.

Usage

1.Obtain and prepare the required forecast data, following the data descriptions in the paper.

2.Configure the input data paths and relevant parameters in the Python script.

3.Run the script to calculate (E_{jk}).

The resulting (E_{jk}) values can be used for subsequent forecast error growth and predictability analyses.

Contact

For questions about the code, please contact the corresponding author of the associated paper.
