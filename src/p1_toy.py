"""Get goal values and dvh from optimal solutions."""
import sys

import numpy as np

import connect
import sample

repo_path = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\'
sys.path.append(repo_path + 'src\\')


def get_dvh(plan, roi_names):
    """Get dvh curves from treatment plan.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObjects
        Current treatment plan.
    roi_names : list
        Regions of interest to include in results.

    Returns
    -------
    dict
        Dose and volumes for given regions of interest.

    """
    dose = plan.TreatmentCourse.TotalDose
    max_dose = max([dose.GetDoseStatistic(RoiName=roi, DoseType='Max')
                    for roi in roi_names])
    dvh = {'Dose': np.linspace(0, max_dose, 100)}
    for roi in roi_names:
        dvh[roi] = dose.GetRelativeVolumeAtDoseValues(RoiName=roi,
                                                      DoseValues=dvh['Dose'])
    return dvh


if __name__ == '__main__':

    # File paths
    funcs_path = repo_path + 'results\\patient1\\funcs_toy.csv'
    goals_path = repo_path + 'results\\patient1\\goals.csv'
    save_path = repo_path + 'results\\patient1\\'

    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Define functions and goals
    funcs, goals, roi_names, pars, results, _ = sample.init_prob(funcs_path,
                                                                 goals_path)

    # Rib Parameters
    rib_pars = {('Rib', 'DoseLevel'): [
        3200, # default
        1400.0, # grid search
        1454.3751069922719, # random search
        1384.7830992256718 # bayesian search
    ]}

    # Calculate plans
    for ii in range(4):
        pars = sample.grid_pars(ii, funcs, pars, rib_pars, 4)
        pars.to_pickle(save_path + 'pars_opt.npy')
        sample.set_pars(plan, pars)
        flag = sample.calc_plan(plan, beam_set, 'PTV', 4800, 95)
        results = sample.get_results(plan, ii, flag, goals, results)
        results.to_pickle(save_path + 'results_opt.npy')
        dvh = get_dvh(plan, roi_names)
        np.save(save_path + f'dvh_{ii}.npy', dvh)
