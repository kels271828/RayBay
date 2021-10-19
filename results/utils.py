import numpy as np

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

plan_types = ['clinical', 'default', 'random', 'bayes']

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
    raise Exception(f"Invalid plan_type '{plan_type}'.")


def get_plan_time(patient, plan_type, stop=False):
    """Get planning time."""
    if plan_type in ['clinical', 'default']:
        raise Exception(f"Time not recorded for {plan_type} plans.")
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
    """Get planning time used for `ii` iterations."""
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


def get_log_path(plan_type):
    log_path = get_plan_path(plan_type).replace('res', 'log')
    return log_path.replace('pkl', 'txt')


def get_best_idx(patient, plan_type, stop=False, n=20, m=15, p=1):
    """Get index of best utility with based on stopping condition."""
    if plan_type in ['clinical', 'default']:
        return 0
    plan = np.load(patient + get_plan_path(plan_type), allow_pickle=True)
    util_vec = plan.opt_result.func_vals
    if stop:
        stop_idx = get_stop_idx(util_vec)
        util_vec = util_vec[:stop_idx + 1]
    return np.argmin(util_vec)










