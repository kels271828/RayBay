"""Solve treatment plan with Bayesian optimization."""
import numpy as np
from skopt import gp_minimize
from time import time

import connect


# What if RayStation can't calculate a plan? Need try/except statements.


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
    weight = 0.75
    app_avg = 318.459
    app_max = 6076.125
    oar_dec = (results[0] - app_avg)/app_avg
    ptv_inc = (results[1] - app_max)/app_max
    return (1 - weight)*oar_dec + weight*ptv_inc


def save_results(plan, beam_set, pars, fpath, normalize=True):
    """Save treatment plan results.

    Parameters
    ----------
    plan : connect.connect_cpytyon.PyScriptObject
        Current RayStation treatment plan.
    beam_set : connect.connect_cpytyon.PyScriptObject
        Current RayStation beam set.
    pars : list
        OAR dose level and percent volume parameters.
    fpath : str
        File path to results directory.
    normalize : bool, optional
        If True, normalized dose to PTV D85 = 4800 cGy. 

    """
    plan_dose = plan.PlanOptimizations[0].PlanningPhaseDose

    # Calculate plan
    [oar_avg, ptv_max] = calc_plan(plan, pars, normalize=normalize)
    if normalize:
        beam_set.NormalizeToPrescription(RoiName='PTV', DoseValue=4800.0,
                                         DoseVolume=95.0,
                                         PrescriptionType='DoseAtVolume')
        oar_avg = plan_dose.GetDoseStatistic(RoiName='Lungs',
                                             DoseType='Average')
        ptv_max = plan_dose.GetDoseStatistic(RoiName='PTV', DoseType='Max')
    
    # Get dose-volume histogram curves
    oar_max = plan_dose.GetDoseStatistic(RoiName='Lungs', DoseType='Max')
    max_dose = max(ptv_max, oar_max)
    dvh_dose = np.linspace(0, max_dose, num=100)
    dvh_oar = plan_dose.GetRelativeVolumeAtDoseValues(RoiName='Lungs',
                                                      DoseValues=dvh_dose)
    dvh_ptv = plan_dose.GetRelativeVolumeAtDoseValues(RoiName='PTV',
                                                      DoseValues=dvh_dose)

    # Save results
    fend = ('.npy', '_norm.npy')[normalize]
    np.save(fpath + 'oar_avg' + fend, oar_avg)
    np.save(fpath + 'ptv_max' + fend, ptv_max)
    np.save(fpath + 'dvh_dose' + fend, dvh_dose)
    np.save(fpath + 'dvh_oar' + fend, dvh_oar)
    np.save(fpath + 'dvh_ptv' + fend, dvh_ptv)


if __name__ == '__main__':
    start = time()
    fpath = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\results\\'

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
   #np.save(fpath + 'results.npy', results) # doesn't work!
    
    # Save plan results
    print('\nSaving results')
    save_results(plan, beam_set, results.x, fpath, normalize=False)
    save_results(plan, beam_set, results.x, fpath, normalize=True)
    
    # Print total time
    total_time = time() - start
    print('Total time: {total_time} seconds')
 