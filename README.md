# autoray
Automatic parameter tuning for RayStation

# Current Workflow

Here is my current workflow for running Python scripts for RayStation in Spyder.
Based on [Landon's instructions](https://github.com/kels271828/autoray/blob/master/docs/Raystation%20Virtual%20Environment.docx), but some things have been modified.
Other scripting information can be found in the [RayStation documentation](https://github.com/kels271828/autoray/blob/master/docs/RSL-D-RS-8B-SG-EN-1.0-2018-12-20%20RayStation%208B%20Scripting%20Guideline.pdf).

## Launch RayStation
1. Log into Citrix Receiver at https://access.radonc.washington.edu/
2. Open RayStation 8B SP1
3. Open RayStation Planning
3. Select patient and plan
4. Make sure "Autoscale to Prescription" is turned off
5. Make sure objective function terms are listed in same order as file or data frame
4. In Scripting tab, run script `run_console` (may need to de-select "Show only validated scripts" in Settings)
5. Create a new tab in the RayStation console

## Set up virtual environment

1. Create virtual environment: ``python -m venv fpath\virtual_environment``
2. Activate virtual environment: ``fpath\virtual_environment\Scripts\activate.bat``
3. Install packages: ``pip install matplotlib numpy pandas pythonnet scikit-optimize spyder`` or ``pip install -r requirements.txt``

Note: Setting up the virtual environment only needs to be done once.

## Run Spyder and connect to RayStation
1. Activate virtual environment: ``fpath\virtual_environment\Scripts\activate.bat``
2. Run Spyder: ``spyder3``
3. Start Spyder kernel: ``python -m spyder_kernels.console``
4. Connect Spyder to RayStation:
    * In Spyder, click options icon at the top left of the Spyder console
    * Click "Connect to existing kernel"
    * Browse to the json file created in the RayStation console
    
## Run Python Script
1. Open script in Spyder
    * To use scripts on local computer, may need to change read/write settings first
2. Click "Run file" button, press F5 button, or use ``%run fpath\script.py`` in console

Haven't tried, but assuming could also run in RayStation console with ``python fpath\script``.
Wouldn't be good for debugging, but okay once script working.
