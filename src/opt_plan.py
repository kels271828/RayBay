"""Optimize treatment plans for patient.

In spyder preferences, select default working directory = the current
working directory.

I am getting new error "Qt Qtwebengineprocess has stopped working"
when I open spyder and when I connect to the kernel.

"""
import pickle
import sys

repo = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\'
sys.path.append(repo + 'src\\')
import optimize

# Setup
patient = repo + 'results\\SBRT_lung_minsun\\'
case = 'pars\\'
utility = 'linear_quadratic'
solver = 'gp_minimize'

# Optimize treatment plan
log_path = patient + case + 'log' + '_' + utility + '_' + solver + '.txt'
with open(log_path, 'w') as fp:
    sys.stdout = fp
    result = optimize.get_plan(
        funcs=patient + case + 'funcs.csv',
        norm=('PTV', 4800, 95),
        goals=patient + 'goals.csv',
        utility=utility,
        solver=solver,
        n_calls=100,
        random_state=7,
        n_initial_points=10,
        verbose=True)

# Save results
result_path = patient + case + 'res' + '_' + utility + '_' + solver + '.pkl'
with open(result_path, 'wb') as fp:
    pickle.dump(result, fp)
