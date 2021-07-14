# RayBay

[![DOI](https://zenodo.org/badge/261853893.svg)](https://zenodo.org/badge/latestdoi/261853893)

Automatic parameter tuning for the RayStation treatment planning system

## Code

* [raybay](/src/raybay.py): Class and functions for specifying and saving problem instances.
* [optimize](/src/optimize.py): Hyperparameter optimization functions for RayStation treatment planning system.
* [analyze](/src/analyze.py): Plotting functions to visualize treatment plan results.

## Workflow

### Define clinical goals and objective functions

* Clinical goals are specified in a CSV file with columns Roi, Type, GoalCriteria, AcceptanceLevel, ParameterValue, Weight, and Shape. Valid types include AverageDose, MinDose, MaxDose, MinDvh, and MaxDvh. Valid shapes include linear and linear_quadratic.
* Constituent functions are specified in a CSV file with columns Roi, FunctionType, DoseLevel, PercentVolume, EudParameterA, and Weight. The row index within the table should correspond to the constituent function in the RayStation objective. Fixed parameters should be a single value, tunable parameters should be a list containing the minimum and maximum values, and irrelevant parameters can be left blank.

See [results](/results) for examples.

Once the virtual environment is set up (see below), you can create clinical goal and objective function files with the following steps:
1. Run [get_volume](/src/get_volume.py) script from Spyder to create initial goals.csv file. This will include the names and volumes for all ROIs in the current case. Delete rows that you don't want to include, and fill out the remaining columns based on your clinical goals. Use columns `Volume (cm^3)` and `RoiVolume (cm^3)` to compute the column `Volume (%)`.
2. Run [setup_patient](/src/setup_patient.py) script from local computer to create goals.csv and funcs.csv files formatted for the approved, default, and bayes cases. For the approved plan, the file approved/funcs.csv will be filled in when you run the script [get_plan](/src/get_plan.py). For the optimized plan, the file bayes/goals.csv will have default weights = 1 for all ROIs and shape = 'linear_quadratic' for all ROIs except for the chestwall and ribs. Update `Weight` and `Shape` columns if needed. The file bayes/funcs.csv will have `DoseLevel` = [gamma/4, gamma] for all goals except for the PTV. The PTV D95 will be 4800, and the PTV Max will be [(gamma - 3*4800)/4, gamma]. Update `DoseLevel`, `PercentVolume`, and `Weight` columns if needed (parameter will be tuned if a range of values is given, and parameter will be held constant if a single value is given).  

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

See [src](/src) for example scripts:
* [get_plan](/src/get_plan.py) gets clinical goal values for the approved plan and the default plan. Make sure you have the correct plan open in RayStation. For the default plan, the constituent functions in RayStation will be cleared and re-populated based on default/funcs.csv.
* [opt_plan](/src/opt_plan.py) gets the optimized treatment plan using either random sampling (`dummy_minimize`) or Bayesian optimization (`gp_minimize`).

## Paper

Data and notebooks used to create the figures appearing in our [preprint](https://arxiv.org/abs/2105.07024) "A feasibility study of a hyperparameter tuning approach to automated inverse planning in radiotherapy" can be found in [results](/results).

#### Abstract
Radiotherapy inverse planning requires treatment planners to modify multiple parameters in the objective function to produce clinically acceptable plans. Due to manual steps in this process, plan quality can vary widely depending on planning time available and planner's skills. The purpose of this study is to automate the inverse planning process to reduce active planning time while maintaining plan quality. We propose a hyperparameter tuning approach for automated inverse planning, where a treatment plan utility is maximized with respect to the limit dose parameters and weights of each organ-at-risk (OAR) objective. Using 6 patient cases, we investigated the impact of the choice of dose parameters, random and Bayesian search methods, and utility function form on planning time and plan quality. For given parameters, the plan was optimized in RayStation, using the scripting interface to obtain the dose distributions deliverable. We normalized all plans to have the same target coverage and compared the OAR dose metrics in the automatically generated plans with those in the manually generated clinical plans. Using 100 samples was found to produce satisfactory plan quality, and the average planning time was 2.3 hours. The OAR doses in the automatically generated plans were lower than the clinical plans by up to 76.8%. When the OAR doses were larger than the clinical plans, they were still between 0.57% above and 98.9% below the limit doses, indicating they are clinically acceptable. For a challenging case, a dimensionality reduction strategy produced a 92.9% higher utility using only 38.5% of the time needed to optimize over the original problem. This study demonstrates our hyperparameter tuning framework for automated inverse planning can significantly reduce the treatment planner's planning time with plan quality that is similar to or better than manually generated plans.
