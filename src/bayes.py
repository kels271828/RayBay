"""Solve treatment plan with Bayesian optimization."""
from time import time

import numpy as np
from skopt import gp_minimize

import connect
import raystation


def objective(plan, pars):
    """Objective function for Bayesian optimization.

    Parameters
    ----------
    plan : connect.connect_cpytyon.PyScriptObject
        Current RayStation treatment plan.
    pars : list
        OAR dose level and percent volume parameters.

    Returns
    -------
    float
        Treatment plan score.

    """
    results = calc_plan(plan, pars)
    return score_plan(results)


def calc_plan(plan, pars, normalize=True):
    """Calculate treatment plan.

    Parameters
    ----------
    plan : connect.connect_cpytyon.PyScriptObject
        Current RayStation treatment plan.
    pars : list
        OAR dose level and percent volume parameters.
    normalize : bool, optional
        If True, estimate the dose normalized to PTV D85 = 4800 cGy.

    Returns
    -------
    list
        OAR average and PTV max doses.

    """
    # Get optimization objects
    plan_opt = plan.PlanOptimizations[0]
    plan_dose = plan_opt.PlanningPhaseDose
    oar_term = plan_opt.Objective.ConstituentFunctions[11] # Lungs

    # Calculate treatment plan
    oar_term.DoseFunctionParameters.DoseLevel = pars[0]
    oar_term.DoseFunctionParameters.PercentVolume = pars[1]
    plan_opt.ResetOptimization()
    plan_opt.RunOptimization()

    # Get results
    oar_avg = plan_dose.GetDoseStatistic(RoiName='Lungs', DoseType='Average')
    ptv_max = plan_dose.GetDoseStatistic(RoiName='PTV', DoseType='Max')

    # Estimate normalized results
    if normalize:
        ptv_d95 = plan_dose.GetDoseAtRelativeVolumes(RoiName='PTV',
                                                     RelativeVolumes=[0.95])[0]
        oar_avg *= 4800/ptv_d95
        ptv_max *= 4800/ptv_d95
    print(f'Dose: {pars[0]}, Volume: {pars[1]}, OAR avg: {oar_avg:e}, PTV max: {ptv_max:e}')
    return [oar_avg, ptv_max]


def score_plan(results):
    """Score treatment plan.

    Parameters
    ----------
    results : list
        OAR average and PTV max doses.

    Returns
    -------
    float
        Treatment plan score.

    """
    weight = 0.25
    app_avg = 318.459
    app_max = 6076.125
    oar_dec = (results[0] - app_avg)/app_avg
    ptv_inc = (results[1] - app_max)/app_max
    return (1 - weight)*oar_dec + weight*ptv_inc


def save_results(roi_names, fpath, normalize=False):
    """Save dose statistics and dvh curves from current plan.

    Parameters
    ----------
    roi_names : list
        Regions of interest to include in results.
    fpath : str
        File path to save results.
    normalize : bool, optional
        If True, normalized dose to PTV D95 = 4800 cGy. 

    """
    if normalize:
        beam_set.NormalizeToPrescription(RoiName='PTV', DoseValue=4800.0,
                                         DoseVolume=95.0,
                                         PrescriptionType='DoseAtVolume')
    fend = fend = ('.npy', '_norm.npy')[normalize]
    np.save(fpath + 'stats' + fend, get_stats(roi_names))
    np.save(fpath + 'dvh' + fend, get_dvh(roi_names))


def get_stats(roi_names):
    """Get dose statistics from current plan.

    Parameters
    ----------
    roi_names : list
        Regions of interest to include in results.

    Returns
    -------
    dict
        Dose statistics for given regions of interest.

    """
    stats = {}
    dose = connect.get_current('Plan').TreatmentCourse.TotalDose
    volumes = [0.99, 0.98, 0.95, 0.5, 0.05, 0.02, 0.01]
    volume_names = ['D99', 'D98', 'D95', 'D50', 'D5', 'D2', 'D1']
    for roi in roi_names:
        stats[roi] = {'Min': dose.GetDoseStatistic(RoiName=roi, DoseType='Min'),
                      'Average': dose.GetDoseStatistic(RoiName=roi,
                                                       DoseType='Average'),
                      'Max': dose.GetDoseStatistic(RoiName=roi, DoseType='Max')}
        doses = dose.GetDoseAtRelativeVolumes(RoiName=roi,
                                              RelativeVolumes=volumes)
        for ii in range(len(volumes)):
            stats[roi][volume_names[ii]] = doses[ii]
    return stats


def get_dvh(roi_names):
    """Get dvh curves from current plan.

    Parameters
    ----------
    roi_names : list
        Regions of interest to include in results.

    Returns
    -------
    dict
        Dose and volumes for given regions of interest.

    """
    dose = connect.get_current('Plan').TreatmentCourse.TotalDose
    max_dose = max([dose.GetDoseStatistic(RoiName=roi, DoseType='Max')
                    for roi in roi_names])
    dvh = {'Dose': np.linspace(0, max_dose, 100)}
    for roi in roi_names:
        dvh[roi] = dose.GetRelativeVolumeAtDoseValues(RoiName=roi,
                                                      DoseValues=dvh['Dose'])
    return dvh


if __name__ == '__main__':
    start = time()
    fpath = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\results\\weight_25\\'

    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Get initial point
    print('Getting initial point')
    x0 = [1500.0, 7.5]
    y0 = objective(plan, x0)

    # Optimize
    print('\nStarting optimization')
    results = gp_minimize(lambda pars: objective(plan, pars),
                          dimensions=[(100, 2000), (0, 10)], x0=x0, y0 = y0,
                          random_state=0, n_calls=25)
    np.save(fpath + 'x.npy', results.x)
    np.save(fpath + 'x_iters.npy', results.x_iters)
    np.save(fpath + 'fun.npy', results.fun)
    np.save(fpath + 'func_vals.npy', results.func_vals)
    
    # Save plan results
    print('\nSaving results')
    save_results(plan, beam_set, results.x, fpath, normalize=False)
    save_results(plan, beam_set, results.x, fpath, normalize=True)
    
    # Print total time
    total_time = time() - start
    print(f'Total time: {total_time} seconds')
 