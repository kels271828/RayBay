"""1D Grid search for treatment plans."""
import pickle
import sys

repo = '\\\\client\\E$\\My Drive\\RayBay\\'
sys.path.append(repo + 'src\\')
import optimize

# Patient
patient = repo + 'results\\SBRT_lung_minsun\\'
#patient = repo + 'results\\ZZ_MK_LLungSBRT3778\\'
#patient = repo + 'results\\ZZ_MK_RLungSBRT4076\\'
#patient = repo + 'results\\ZZ_MK_RULungSBRT3796\\'
#patient = repo + 'results\\ZZ_MK_RLSBRT1931\\'

# Case
case = 'weight\\'

# Calculate treatment plans
log_path = patient + case + 'log_weight.txt'
stdout = sys.stdout
with open(log_path, 'w') as fp:
    sys.stdout = fp
    result = optimize.grid_search(
        funcs=patient + case + 'funcs_weight.csv',
        norm=('PTV', 4800, 95),
        goals=patient + case + 'goals_lin.csv',
        n_points=10,
        weight=True)
sys.stdout = stdout

# Save results
result_path = patient + case + 'res_weight.pkl'
with open(result_path, 'wb') as fp:
    pickle.dump(result, fp)
