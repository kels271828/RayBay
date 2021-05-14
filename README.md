# RayBay

[![DOI](https://zenodo.org/badge/261853893.svg)](https://zenodo.org/badge/latestdoi/261853893)

Automatic parameter tuning for the RayStation treatment planning system

## Code

* [raybay](/src/raybay.py): Class and functions for specifying and saving problem instances.
* [optimize](/src/optimize.py): Hyperparameter optimization functions for RayStation treatment planning system.
* [analyze](/src/analyze.py): Plotting functions to visualize treatment plan results.

Data and notebooks used to create the figures appearing in our preprint "A feasibility study of a hyperparameter tuning approach to automated inverse planning in radiotherapy" can be found in [results](/results).

## Workflow

### Define clinical goals and objective functions

* Clinical goals are specified in a CSV file with columns Roi, Type, GoalCriteria, AcceptanceLevel, ParameterValue, Weight, and Shape. Valid types include AverageDose, MinDose, MaxDose, MinDvh, and MaxDvh. Valid shapes include linear and linear_quadratic.
* Constituent functions are specified in a CSV file with columns Roi, FunctionType, DoseLevel, PercentVolume, EudParameterA, and Weight. The row index within the table should correspond to the constituent function in the RayStation objective. Fixed parameters should be a single value, tunable parameters should be a list containing the minimum and maximum values, and irrelevant parameters can be left blank.

See [results](/results) for examples.

### Launch RayStation
1. Log in to Citrix Receiver
2. Open RayStation 8B SP1
3. Open RayStation Planning
3. Select patient and plan
4. Make sure "Autoscale to Prescription" is turned off
5. Make sure objective function terms are listed in same order as csv file
4. In Scripting tab, run script `run_console` (may need to de-select "Show only validated scripts" in Settings)
5. Create a new tab in the RayStation console

### Set up virtual environment

1. Create virtual environment: ``python -m venv fpath\virtual_environment``
2. Activate virtual environment: ``fpath\virtual_environment\Scripts\activate.bat``
3. Install packages: ``pip install -r requirements.txt``

Note: Setting up the virtual environment only needs to be done once.

### Run Spyder and connect to RayStation
1. Activate virtual environment: ``fpath\virtual_environment\Scripts\activate.bat``
2. Run Spyder: ``spyder``
3. Start Spyder kernel: ``python -m spyder_kernels.console``
4. Connect Spyder to RayStation:
    * In Spyder, click options icon at the top left of the Spyder console
    * Click "Connect to existing kernel"
    * Enter kernel number listed in the RayStation console
    
### Run Python Script
1. Open script in Spyder
    * To use scripts on local computer, may need to change read/write settings first
2. Click "Run file" button

See [src](/src) for example scripts.
