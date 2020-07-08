"""Optimize treatment plans for patient."""
import sys

import opt

repo_path = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\'
sys.path.append(repo_path + 'src\\')
funcs_path = repo_path + 'results\\patient1\\funcs_toy.csv'
goals_path = repo_path + 'results\\patient1\\goals.csv'
save_path = repo_path + 'results\\patient1\\'
#opt.gp_minimize(funcs_path, 'PTV', 4800, 95, goals_path, save_path)
#opt.forest_minimize(funcs_path, 'PTV', 4800, 95, goals_path, save_path)
opt.dummy_minimize(funcs_path, 'PTV', 4800, 95, goals_path, save_path)
