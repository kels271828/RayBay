# autoray
Automatic parameter tuning for RayStation

# Current Workflow

Here is my current workflow for running Python scripts for RayStation in Spyder.
Based on [Landon's instructions](https://github.com/kels271828/autoray/blob/master/Raystation%20Virtual%20Environment.docx), but some things have been modified.

## Launch RayStation
1. Log into RayStation
2. Select patient
3. In Scripting tab, run script `run_console` (may need to de-select "Show only validated scripts" in Settings)
4. Create a new tab in the RayStation console

## Set up virtual environment

1. Create virtual environment: ``python -m venv fpath\virtual_environment``
2. Activate virtual environment: ``fpath\virtual_environment\Scripts\activate.bat``
3. Install packages: ``pip install matplotlib numpy pythonnet spyder`` or ``pip install -r requirements.txt``

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
2. Click "Run file" button, press F5 button, or use ``%run 'fpath\script.py`` in console
