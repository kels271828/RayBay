"""Get default and approved plans.

Current workflow:
    * In spyder preferneces, select
        Default working directory: The current working directory

"""
import pickle
import sys

import connect

repo_path = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\'
sys.path.append(repo_path + 'src\\')
import optimize
import raybay


# Setup
patient_path = repo_path + 'results\\ZZ_MK_LLungSBRT3778\\'
case_path = 'default\\'

# Get RayStation objects
patient = connect.get_current('Patient')
case = connect.get_current('Case')
plan = connect.get_current('Plan')
beam_set = connect.get_current('BeamSet')

# Initialize result object
result = raybay.RaybayResult(
    patient.Name,
    case.CaseName,
    plan.Name,
    patient_path + case_path + 'funcs.csv',
    ('PTV', 4800, 95),
    'N/A',
    patient_path + 'goals.csv')

# Add results
if 'default' in case_path:
    optimize.set_pars(plan, result.funcs, [])
    flag = optimize.calc_plan(plan, beam_set, result.norm)
else:
    flag = 0
result.opt_result = optimize.get_score(plan, result.goals, flag, 
                                       result.goal_result)
result.dvh_result = optimize.get_dvh(result.roi_list)

# Save results
with open(patient_path + case_path + 'result', 'wb') as fp:
    pickle.dump(result, fp)
