"""Get clinical and default plans."""
import pickle
import sys

import connect

repo_path = '\\\\client\\E$\\My Drive\\RayBay\\'
sys.path.append(repo_path + 'src\\')
import optimize
import raybay

# Patient
patient_path = repo_path + 'results\\SBRT_lung_minsun\\'
#patient_path = repo_path + 'results\\ZZ_MK_LLungSBRT3778\\'
#patient_path = repo_path + 'results\\ZZ_MK_RLungSBRT4076\\'
#patient_path = repo_path + 'results\\ZZ_MK_RULungSBRT3796\\'
#patient_path = repo_path + 'results\\ZZ_MK_RLSBRT1931\\'
#patient_path = repo_path + 'results\\ZZ_MK_LLLungSBRT2736\\'
#patient_path = repo_path + 'results\\ZZ_MK_LULSBRT4544\\'
#patient_path = repo_path + 'results\\ZZ_MK_SBRTLL0924allviolated\\'
#patient_path = repo_path + 'results\\ZZ_MK_SBRTLL7289\\'
#patient_path = repo_path + 'results\ZZ_MK_SBRTLLL8973\\'
#patient_path = repo_path + 'results\\ZZ_MK_SBRTRL7289\\'
#patient_path = repo_path + 'results\\ZZ_MK_SBRTRUL_2928allviolate\\'

# Case
case_path = 'approved\\'
#case_path = 'default\\'

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
    ('PTV_4800', 4800, 95),  # check PTV name
    goals=patient_path + case_path + 'goals.csv')

# Add results
if 'default' in case_path:
    optimize.add_funcs(plan, result.func_df)
    optimize.set_pars(plan, result.func_df, [])
    flag = optimize.calc_plan(plan, beam_set, result.norm)
else:
    result.func_df = optimize.get_funcs(plan)
    result.func_df.to_csv(patient_path + case_path + 'funcs.csv')
    flag = 0
result.flag_list.append(flag)
goal_results = optimize.get_results(plan, result.goal_df)
scale = optimize.get_scale(result.goal_df, result.norm, goal_results)
for index, row in result.goal_df.iterrows():
    result.goal_dict[index].append(scale*goal_results[index])
result.dvh_dict = optimize.get_dvh(result.roi_list)
coeff = result.goal_dict[6][0]/result.dvh_dict['Dose'][-1]
result.dvh_dict['Dose'] *= coeff  # normalize (check index in prev line)

# Save results
result_path = patient_path + case_path + 'res_' + case_path[:-1] + '.pkl'
with open(result_path, 'wb') as fp:
    pickle.dump(result, fp)
