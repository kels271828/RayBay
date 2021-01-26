"""1D Grid search for treatment plans.

In spyder preferences, select default working directory = the current
working directory.

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
case = 'grid\\'

# Calculate treatment plans


# Calculate treatment plans
log_path = patient + case + 'log_grid_1.txt'
stdout = sys.stdout
with open(log_path, 'w') as fp:
    sys.stdout = fp
    result = optimize.grid_search(
        funcs=patient + case + 'funcs.csv',
        norm=('PTV', 4800, 95),
        goals=patient + case + 'goals_lin_rib.csv',
        n_points=25)
sys.stdout = stdout

# Save results
result_path = patient + case + 'res_grid_1.pkl'
with open(result_path, 'wb') as fp:
    pickle.dump(result, fp)
