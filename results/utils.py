import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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


par_names = goal_names[:5] + goal_names[6:]


def get_plan(patient, plan_type):
    return np.load(patient + get_plan_path(plan_type), allow_pickle=True)


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
    return 100*(row[value] - row[reference])/np.abs(row[reference])


### Time Results ###


def get_time_df(plan_type, stop=False):
    """Get planning time for all patients."""
    df = pd.DataFrame({
        'patient': patients,
        'plan_type': len(patients)*[plan_type],
        'plan_time': [get_plan_time(patient, plan_type, stop)
                      for patient in patients],
        'plan_iter': [get_plan_iter(patient, plan_type, stop)
                      for patient in patients]})
    return df


def get_plan_time(patient, plan_type, stop=False):
    """Get planning time."""
    plan = get_plan(patient, plan_type)
    if stop:
        util_vec = plan.opt_result.func_vals
        ii = get_stop_idx(util_vec)
        return get_log_time(ii, patient, plan_type)
    return plan.time/3600.0


def get_plan_iter(patient, plan_type, stop=False):
    if stop:
        plan = get_plan(patient, plan_type)
        util_vec = plan.opt_result.func_vals
        ii = get_stop_idx(util_vec)
        return ii + 1
    return 100


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
    """Get utility for all patients."""
    df = pd.DataFrame({
        'patient': patients,
        'plan_type': len(patients)*[plan_type],
        'plan_util': [get_plan_util(patient, plan_type, stop)
                      for patient in patients]})
    return df


def get_plan_util(patient, plan_type, stop=False):
    """Get plan utility."""
    plan = get_plan(patient, plan_type)
    if plan_type in ['clinical', 'default']:
        ref_plan = np.load(patient + get_plan_path('random'), allow_pickle=True)
        return raybay.get_utility(ref_plan.goal_df, plan.goal_dict)[0]
    if stop:
        ii = get_best_idx(patient, plan_type, stop=True)
        util_vec = plan.opt_result.func_vals
        return -util_vec[ii]
    return -plan.opt_result.fun


def get_best_idx(patient, plan_type, stop=False, n=20, m=15, p=1):
    """Get index of best utility based on stopping condition."""
    if plan_type in ['clinical', 'default']:
        return 0
    plan = get_plan(patient, plan_type)
    util_vec = plan.opt_result.func_vals
    if stop:
        stop_idx = get_stop_idx(util_vec)
        util_vec = util_vec[:stop_idx + 1]
    return np.argmin(util_vec)


### Parameter Results ###


def get_pars_df(plan_type, stop=False):
    """Get plan parameters for all patients."""
    df = pd.concat([get_plan_pars(patient, plan_type, stop)
                    for patient in patients])
    return df


def get_plan_pars(patient, plan_type, stop=False):
    """Get plan parameters."""
    goal_vals = get_goal_vals(patient, plan_type)
    df = pd.DataFrame({
        'patient': len(par_names)*[patient],
        'plan_type': len(par_names)*[plan_type],
        'par_name': par_names,
        'par_val': get_par_vals(patient, plan_type, stop),
        'goal_val': goal_vals[:5] + goal_vals[6:]})
    return df


def get_par_vals(patient, plan_type, stop=False):
    """Get vector of plan parameters."""
    plan = get_plan(patient, plan_type)
    if stop:
        ii = get_best_idx(patient, plan_type, stop=True)
        x_iters = plan.opt_result.x_iters
        return x_iters[ii]
    return plan.opt_result.x


def get_goal_vals(patient, plan_type):
    """Get plan clinical goal values."""
    plan = get_plan(patient, plan_type)
    return plan.goal_df['AcceptanceLevel'].tolist()


### Dose Results ###


def get_dose_df(plan_type, stop=False):
    """Get plan dose values for all patients."""
    df = pd.concat([get_plan_dose(patient, plan_type, stop)
                    for patient in patients])
    return df


def get_plan_dose(patient, plan_type, stop=False):
    """Get plan dose values."""
    goal_vals = get_goal_vals(patient, plan_type)
    df = pd.DataFrame({
        'patient': len(par_names)*[patient],
        'plan_type': len(par_names)*[plan_type],
        'dose_name': par_names,
        'dose_val': get_dose_vals(patient, plan_type, stop),
        'goal_val': goal_vals[:5] + goal_vals[6:]})
    return df


def get_dose_vals(patient, plan_type, stop=False):
    """Get vector of plan dose values."""
    plan = get_plan(patient, plan_type)
    ii = get_best_idx(patient, plan_type, stop)
    dose_vals = [plan.goal_dict[goal][ii] for goal in plan.goal_dict]
    return dose_vals[:5] + dose_vals[6:]


def heatmap(df, col, col_types, diff_type, label):
    for patient in patients:
        df_sub = df[df['patient'] == patient]
        dose_vals = np.array([df_sub[df_sub[col] == col_type][diff_type].values
                              for col_type in col_types])
        fig, ax = plt.subplots(figsize=(dose_vals.shape[1], dose_vals.shape[0]))
        sns.heatmap(dose_vals, cmap=sns.diverging_palette(220, 20, n=256),
                    center=0, annot=True, fmt=".2f", ax=ax,
                    cbar_kws={'label': f"Percent Difference from {label}"})
        ax.set_xticklabels(par_names, rotation=90)
        ax.set_yticklabels(col_types, rotation=0)
        ax.set_title(patient)
