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
