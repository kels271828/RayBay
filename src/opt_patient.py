"""Optimize treatment plans for patient."""
import pickle
import sys

import optimize

repo_path = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\'
sys.path.append(repo_path + 'src\\')
funcs_path = repo_path + 'results\\SBRT_lung_minsun\\funcs_small.csv'
goals_path = repo_path + 'results\\SBRT_lung_minsun\\goals.csv'
save_path = repo_path + 'results\\SBRT_lung_minsun\\small\\'
result = optimize.plan_opt(funcs_path, ('PTV', 4800, 95), goals_path,
                           n_calls=5, n_initial_points=2)
pickle.dump(result, open(save_path + 'results', 'w'))
