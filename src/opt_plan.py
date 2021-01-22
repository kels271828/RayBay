"""Optimize treatment plans for patient.

In spyder preferences, select default working directory = the current
working directory.

I am getting new error "Qt Qtwebengineprocess has stopped working"
when I open spyder and when I connect to the kernel.

Be careful about the goals csv now, since I've added columns for utility!

"""
import pickle
import sys

repo = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\'
sys.path.append(repo + 'src\\')
import optimize

# Patient
patient = repo + 'results\\SBRT_lung_minsun\\'
#patient = repo + 'results\\ZZ_MK_LLungSBRT3778\\'
#patient = repo + 'results\\ZZ_MK_RLungSBRT4076\\'
#patient = repo + 'results\\ZZ_MK_RULungSBRT3796\\'

# Case
#case = 'sample\\'
#case = 'dummy\\'
case = 'bayes\\'

# Solver
#solver = 'dummy_minimize'
solver = 'gp_minimize'
#solver = 'forest_minimize'

# Optimize treatment plan
log_path = patient + case + 'log_' + solver + '.txt'
stdout = sys.stdout
with open(log_path, 'w') as fp:
    sys.stdout = fp
    result = optimize.get_plan(
        funcs=patient + case + 'funcs.csv',
        norm=('PTV', 4800, 95),
        goals=patient + case + 'goals.csv',
        solver=solver,
        n_calls=50,
        random_state=1,
        n_initial_points=10,
        verbose=True)
sys.stdout = stdout

# Save results
result_path = patient + case + 'res_' + solver + '.pkl'
with open(result_path, 'wb') as fp:
    pickle.dump(result, fp)
