"""Sample treatment plans for patient 2."""
import sys

import sample

repo_path = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\'
sys.path.append(repo_path + 'src\\')
funcs_path = repo_path + 'results\\patient2\\funcs_half.csv'
goals_path = repo_path + 'results\\patient2\\goals.csv'
save_path = repo_path + 'results\\patient2\\'
sample.sample_plans(funcs_path, 'PTV 4/7/20', 6270, 99, goals_path, repo_path)
