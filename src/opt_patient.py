"""Optimize treatment plans for patient."""
import dill
import sys

import optimize

repo_path = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\'
sys.path.append(repo_path + 'src\\')
funcs_path = repo_path + 'results\\SBRT_lung_minsun\\funcs_small.csv'
goals_path = repo_path + 'results\\SBRT_lung_minsun\\goals.csv'
save_path = repo_path + 'results\\SBRT_lung_minsun\\'
result = optimize.plan_opt(funcs_path, ('PTV', 4800, 95), goals_path,
                           n_calls=3, n_initial_points=2)
dill.dump(result, open(save_path + 'result', 'wb'))

# can't save local objects... nested functions...
# result.opt_result.specs['args']['funcs'] = 'local'
# removes what was 
# <function optimize.plan_opt.<locals>.obj(pars)>
# which was preventing it from being saved...
# 
# would pickle also still work?
# yes, don't need dill anymore...

# can't figure out how else to pass function to gp_minimize
# other than inner function or lambda...