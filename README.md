# autoray
Automatic parameter tuning for RayStation

# Connecting Spyder to RayStation

## Launch RayStation
1. Log into RayStation
2. Select patient
3. In Scripting tab, run script `run_console` (may need to de-select "Show only validated scripts" in Settings)

## Set up virtual environment
1. Create a new tab in the RayStation console
2. Create virtual environment: ``python -m venv C:\Users\kmaass\Documents\virtual_environment
3. Activate virtual environment: ``C:\Users\kmaass\Documents\virtual_environment\Scripts\activate.bat``
4. Install packages: ``pip install numpy spyder pythonnet``

## Run Spyder and connect to RayStation
1. Start Spyder kernel: ``python -m spyder_kernels.console``
2. Create a new tab in the RayStation console
3. Activate virtual environment: ``C:\Users\kmaass\Documents\virtual_environment\Scripts\activate.bat``
4. Run Spyder: ``spyder3``
5. Connect to RayStation:
    * Click gear icon at the top left of the Spyder console
    * Click "Connect to existing kernel"
    * Browse to the json file created in the RayStation console
    
# Notes
    
At the moment, the virtual environments I've created aren't saved after I log out of RayStation.
One thing I might try is to make a bash script I can copy and paste into the console to automate the above steps.
How do I save the plots or results I generate from a script? Can I use scp?
