"""Optimize treatment plans for patient.

In spyder preferences, select default working directory = the current
working directory.

I am getting new error "Qt Qtwebengineprocess has stopped working"
when I open spyder and when I connect to the kernel.

"""
import pickle
import sys

import optimize

repo = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\'
sys.path.append(repo + 'src\\')

# Setup (CHECK RAYSTATION PATIENT, CASE, AND SOLVER!)
patient = repo + 'results\\SBRT_lung_minsun\\'
case = 'dummy\\'

# Optimize treatment plan
with open(patient + case + 'log.txt', 'w') as fp:
    sys.stdout = fp
    result = optimize.get_plan(
        funcs=patient + case + 'funcs2.csv',
        norm=('PTV', 4800, 95),
        goals=patient + 'goals.csv',
        solver='dummy_minimize',
        n_calls=50,
        random_state=7,
        n_initial_points=10,
        verbose=True)

# Save results
with open(patient + case + 'result2', 'wb') as fp:
    pickle.dump(result, fp)
