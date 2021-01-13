"""Get approved and default plans.

In Spyder preferences, set current working directory to default working
directory.

"""
import pickle
import sys

import connect

repo_path = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\'
sys.path.append(repo_path + 'src\\')
import optimize
import raybay

# Patient
patient_path = repo_path + 'results\\SBRT_lung_minsun\\'
#patient_path = repo_path + 'results\\ZZ_MK_LLungSBRT3778\\'
#patient_path = repo_path + 'results\\ZZ_MK_RLungSBRT4076\\'
#patient_path = repo_path + 'results\\ZZ_MK_RULungSBRT3796\\'

# Case
case_path = 'approved\\'
#case_path = 'default\\'

# Get RayStation objects
patient = connect.get_current('Patient')
case = connect.get_current('Case')
plan = connect.get_current('Plan')
beam_set = connect.get_current('BeamSet')

# Initialize result object
result = raybay.RaybayResult(patient.Name, case.CaseName, plan.Name,
                             patient_path + case_path + 'funcs.csv',
                             ('PTV', 4800, 95), patient_path + 'goals.csv')

# Add results
if 'default' in case_path:
    optimize.set_pars(plan, result.func_df, [])
    flag = optimize.calc_plan(plan, beam_set, result.norm)
else:
    flag = 0
result.flag_list.append(flag)
goal_results = optimize.get_results(plan, result.goal_df)
for index, row in result.goal_df.iterrows():
    result.goal_dict[index].append(goal_results[index])
result.dvh_dict = optimize.get_dvh(result.roi_list)

# Save results
with open(patient_path + case_path + 'result', 'wb') as fp:
    pickle.dump(result, fp)
