"""Optimize treatment plans."""
import pickle
import sys

repo = '\\\\client\\E$\\My Drive\\RayBay\\'
sys.path.append(repo + 'src\\')
import optimize

# Patient
#patient = repo + 'results\\SBRT_lung_minsun\\'
#patient = repo + 'results\\ZZ_MK_LLungSBRT3778\\'
#patient = repo + 'results\\ZZ_MK_RLungSBRT4076\\'
#patient = repo + 'results\\ZZ_MK_RULungSBRT3796\\'
#patient = repo + 'results\\ZZ_MK_RLSBRT1931\\'
#patient = repo + 'results\\ZZ_MK_LLLungSBRT2736\\'
patient = repo + 'results\\ZZ_MK_LULSBRT4544\\'

# Case
case = 'bayes\\'
#case = 'dimension\\'

# Solver
#solver = 'dummy_minimize'
solver = 'gp_minimize'

# Optimize treatment plan
log_path = patient + case + 'log_linquad_' + solver + '.txt'
stdout = sys.stdout
with open(log_path, 'w') as fp:
    sys.stdout = fp
    result = optimize.get_plan(
        funcs=patient + case + 'funcs.csv',
        norm=('PTV', 4800, 95),
        goals=patient + case + 'goals_linquad.csv',
        solver=solver,
        n_calls=100,
        random_state=1,
        n_initial_points=20,
        verbose=True)
sys.stdout = stdout

# Save results
result_path = patient + case + 'res_linquad_' + solver + '.pkl'
with open(result_path, 'wb') as fp:
    pickle.dump(result, fp)
