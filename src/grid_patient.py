"""Compute grid search of treatment plans for patient."""
import sys

import sample

repo_path = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\'
sys.path.append(repo_path + 'src\\')
funcs_path = repo_path + 'results\\patient1\\funcs_grid.csv'
goals_path = repo_path + 'results\\patient1\\goals.csv'
save_path = repo_path + 'results\\patient1\\'
sample.grid_search(funcs_path, 'PTV', 4800, 95, goals_path, save_path,
                   n_points=25)
