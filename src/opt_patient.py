"""Optimize treatment plans for patient.

Current workflow:
    * In spyder preferences, select
        Default working directory: The directory of the file being executed
    * Import `connect` module while still in working directory
        C:\Program Files\RaySearch Laboratories\RayStation 8B SP1\ScriptClient
      before running the script.

I am getting new error
    "Qt Qtwebengineprocess has stopped working"
when I open spyder and when I connect to the kernel.
Not sure if it's important?

"""
import pickle
import sys

import optimize

# Setup
patient = '..\\results\\SBRT_lung_minsun\\'
case = 'bayes_pars_2\\'

# Optimize treatment plan
with open(patient + case + 'log.txt', 'w') as fp:
    sys.stdout = fp
    result = optimize.get_plan(
        funcs=patient + case + 'funcs.csv',
        norm=('PTV', 4800, 95),
        goals=patient + 'goals.csv',
        solver='gp_minimize',
        n_calls=50,
        random_state=7,
        n_initial_points=10,
        verbose=True
    )

# Save results
with open(patient + case + 'result', 'wb') as fp:
    pickle.dump(result, fp)
