import sys

import numpy as np
import pandas as pd

sys.path.append('../src')
import raybay


patients = [
    'SBRT_lung_minsun',
    'ZZ_MK_LLungSBRT3778',
    'ZZ_MK_LULSBRT4544',
    'ZZ_MK_RLSBRT1931',
    'ZZ_MK_RLungSBRT4076',
    'ZZ_MK_RULungSBRT3796',
    'ZZ_MK_SBRTLL7289',
    'ZZ_MK_SBRTLLL8973',
    'ZZ_MK_SBRTRL7289',
    'ZZ_MK_SBRTRUL_2928allviolate'
]


goal_names = [
    '1_Chestwall_MaxDVH',
    '2_D2cm_MaxDose',
    '3_Esophagus_MaxDVH',
    '4_Lungs_MaxDVH',
    '5_Lungs_MaxDVH',
    '6_PTV_MinDVH',
    '7_PTV_MaxDose',
    '8_Rib_MaxDVH',
    '9_Rib_MaxDose',
    '10_SpinalCord_MaxDVH',
    '11_SpinalCord_MaxDose'
]


def get_plan_path(plan_type):
    if plan_type == 'clinical':
        return '/approved/res_approved.pkl'
    if plan_type == 'default':
        return '/default/res_default.pkl'
    if plan_type == 'random':
        return '/bayes/res_linquad_dummy_minimize.pkl'
    if plan_type == 'bayes':
        return '/bayes/res_linquad_gp_minimize.pkl'


def get_log_path(plan_type):
    log_path = get_plan_path(plan_type).replace('res', 'log')
    return log_path.replace('pkl', 'txt')


def get_percent_diff(row, value, reference):
    return 100*(row[value] - row[reference])/row[reference]


### Time Results ###


def get_time_df(plan_type, stop=False):
    df = pd.DataFrame({
        'patient': patients,
        'plan_type': len(patients)*[plan_type],
        'plan_time': [get_plan_time(patient, plan_type, stop)
                      for patient in patients]})
    return df


def get_plan_time(patient, plan_type, stop=False):
    """Get planning time."""
    if stop:
        return get_stop_time(patient, plan_type)
    plan = np.load(patient + get_plan_path(plan_type), allow_pickle=True)
    return plan.time/3600.0


def get_stop_time(patient, plan_type):
    """Get planning time with stopping condition."""
    plan = np.load(patient + get_plan_path(plan_type), allow_pickle=True)
    util_vec = plan.opt_result.func_vals
    ii = get_stop_idx(util_vec)
    return get_log_time(ii, patient, plan_type)


def get_stop_idx(util_vec, n=20, m=15, p=1):
    """Get index of last iteration based on stopping condition."""
    best_util = np.minimum.accumulate(util_vec)
    for ii in range(n - 1, len(best_util)):
        max_util = best_util[ii - m + 1]
        min_util = best_util[ii]
        if 100*np.abs((min_util - max_util)/max_util) < p:
            return ii
    return len(best_util) - 1


def get_log_time(ii, patient, plan_type):
    """Get planning time for `ii` iterations."""
    with open(patient + get_log_path(plan_type)) as f:
        log = f.readlines()
    count = 0
    total_time = 0
    for row in log:
        if 'Time' in row:
            total_time += float(row.split()[-1])
            count += 1
        if count > ii:
            return total_time/3600.0
    return total_time/3600.0


### Utility Results ###


def get_util_df(plan_type, stop=False):
    df = pd.DataFrame({
        'patient': patients,
        'plan_type': len(patients)*[plan_type],
        'plan_util': [get_plan_util(patient, plan_type, stop)
                      for patient in patients]})
    return df

def get_plan_util(patient, plan_type, stop=False):
    """Get plan utility."""
    if plan_type in ['clinical', 'default']:
        plan = np.load(patient + get_plan_path(plan_type), allow_pickle=True)
        ref_plan = np.load(patient + get_plan_path('random'), allow_pickle=True)
        return raybay.get_utility(ref_plan.goal_df, plan.goal_dict)[0]
    if stop:
        return get_stop_util(patient, plan_type)
    plan = np.load(patient + get_plan_path(plan_type), allow_pickle=True)
    return -plan.opt_result.fun


def get_stop_util(patient, plan_type):
    """Get plan utility with stopping condition."""
    ii = get_best_idx(patient, plan_type, stop=True)
    plan = np.load(patient + get_plan_path(plan_type), allow_pickle=True)
    util_vec = plan.opt_result.func_vals
    return -util_vec[ii]


def get_best_idx(patient, plan_type, stop=False, n=20, m=15, p=1):
    """Get index of best utility based on stopping condition."""
    if plan_type in ['clinical', 'default']:
        return 0
    plan = np.load(patient + get_plan_path(plan_type), allow_pickle=True)
    util_vec = plan.opt_result.func_vals
    if stop:
        stop_idx = get_stop_idx(util_vec)
        util_vec = util_vec[:stop_idx + 1]
    return np.argmin(util_vec)
