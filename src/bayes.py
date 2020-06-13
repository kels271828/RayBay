"""Solve treatment plan with Bayesian optimization."""
from time import time

import numpy as np
from skopt import gp_minimize

import connect

def objective(plan, beam_set, roi_names, pars):
    """Objective function for Bayesian optimization.

    Parameters
    ----------
    plan : connect.connect_cpytyon.PyScriptObject
        Current treatment plan.
    beam_set : connect.connect_cpytyon.PyScriptObject
        Current beam set.
    roi_names : list
        Regions of interest to include in plan.
    pars : list
        OAR dose level and percent volume parameters.

    Returns
    -------
    float
        Treatment plan score.

    """
    set_pars(plan, pars)
    norm = calc_plan(plan, beam_set)
    if norm:
        return score_plan(plan, roi_names)
    else:
        return 1e2


def set_pars(plan, pars):
    """Set objective function parameters.

    Hardcoded for now:
    pars = [
        1, PTV, MinDose, DoseLevel, 6200
        2, Rib, MaxDvh, DoseLevel, 3200
        3, SpinalCanal, MaxDvh, DoseLevel, 2080
        4, Heart, MaxDvh, DoseLevel, 2800
        5, Chestwall_L, MaxDvh, DoseLevel, 3000
        6, Lungs, MaxDvh, DoseLevel, 2000
    ]

    Parameters
    ----------
    plan : connect.connect_cpytyon.PyScriptObject
        Current treatment plan.
    pars : list
        Objective function parameters.

    """
    funcs = plan.PlanOptimizations[0].Objective.ConstituentFunctions
    roi_idx = [1, 2, 3, 4, 5, 6]
    for ii in range(len(roi_idx)):
        funcs[roi_idx[ii]].DoseFunctionParameters.DoseLevel = pars[ii]


def calc_plan(plan, beam_set):
    """Calculate and normalize treatment plan.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    beam_set : connect.connect_cpython.PyScriptObject
        Current beam set.
        
    Results
    -------
    bool
        Normalization status.

    """
    # Calculate plan
    plan.PlanOptimizations[0].ResetOptimization()
    try:
        plan.PlanOptimizations[0].RunOptimization()
    except:
        return False

    # Normalize plan
    try:
        beam_set.NormalizeToPrescription(RoiName='PTV', DoseValue=4800.0,
                                         DoseVolume=95.0,
                                         PrescriptionType='DoseAtVolume')
        return True
    except:
        return False
        

def score_plan(plan, roi_names):
    """Score treatment plan.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    roi_names : list
        Regions of interest to include in results.

    Returns
    -------
    float
        Treatment plan score.
        
        
        assuming that reason GetClinicalGoalValue() wouldn't work
        is that calc_dose didn't work, so already 1e6

    """
    score = 0
    vals = get_objs(plan)
    for roi in roi_names:
        level = vals[roi][-1]['DoseValue']
        value = vals[roi][-1]['ResultValue']
        score += (value - level)/level
    return score


def get_objs(plan):
    """Get values related to objective terms.

    Totally hard-coded for now.

    Parameters
    ----------
    dose : connect.connect_cpython.PyScriptObject
        Total dose for current treatment plan.

    Returns
    -------
    dict
        Values related to objective terms.

    """
    dose = plan.TreatmentCourse.TotalDose
    values = {}
    values['PTV'] = [{
        'FunctionType': 'MinDose',
        'DoseValue': 4800,
        'PercentVolume': 0,
        'ResultValue': dose.GetDoseStatistic(RoiName='PTV', DoseType='Min')
    }, {
        'FunctionType': 'MaxDose',
        'DoseValue': 6200,
        'PercentVolume': 0,
        'ResultValue': dose.GetDoseStatistic(RoiName='PTV', DoseType='Max')
    }]
    values['Rib'] = [{
        'FunctionType': 'MaxDvh',
        'DoseValue': 3200, 
        'PercentVolume': 0.27,
        'ResultValue': dose.GetDoseAtRelativeVolumes(RoiName='Rib',
            RelativeVolumes=[0.0027])[0]
    }]
    values['SpinalCanal'] = [{
        'FunctionType': 'MaxDvh',
        'DoseValue': 2080,
        'PercentVolume': 0.67,
        'ResultValue': dose.GetDoseAtRelativeVolumes(RoiName='SpinalCanal',
            RelativeVolumes=[0.0067])[0]
    }]
    values['Heart'] = [{
        'FunctionType': 'MaxDvh',
        'DoseValue': 2800,
        'PercentVolume': 1.84,
        'ResultValue': dose.GetDoseAtRelativeVolumes(RoiName='Heart',
            RelativeVolumes=[0.0184])[0]
    }]
    values['Chestwall_L'] = [{
        'FunctionType': 'MaxDvh',
        'DoseValue': 3000,
        'PercentVolume': 2.04,
        'ResultValue': dose.GetDoseAtRelativeVolumes(RoiName='Chestwall_L',
            RelativeVolumes=[0.0204])[0]
    }]
    values['Lungs'] = [{
        'FunctionType': 'MaxDvh',
        'DoseValue': 2000,
        'PercentVolume': 10,
        'ResultValue': dose.GetDoseAtRelativeVolumes(RoiName='Lungs',
            RelativeVolumes=[0.1])[0]
    }]
    return values


def save_results(plan, roi_names, fpath):
    """Save dose statistics and dvh curves from treatment plan.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    roi_names : list
        Regions of interest to include in results.
    fpath : str
        File path to save results.

    """
    np.save(fpath + 'stats.npy', get_stats(plan, roi_names))
    np.save(fpath + 'dvh.npy', get_dvh(plan, roi_names))


def get_stats(plan, roi_names):
    """Get dose statistics from treatment plan.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    roi_names : list
        Regions of interest to include in results.

    Returns
    -------
    dict
        Dose statistics for given regions of interest.

    """
    stats = {}
    dose = plan.TreatmentCourse.TotalDose
    volumes = [0.99, 0.98, 0.95, 0.9, 0.5, 0.1, 0.05, 0.02, 0.01]
    volume_names = ['D99', 'D98', 'D95', 'D90', 'D50', 'D10', 'D5', 'D2', 'D1']
    for roi in roi_names:
        stats[roi] = {'Min': dose.GetDoseStatistic(RoiName=roi, DoseType='Min'),
                      'Max': dose.GetDoseStatistic(RoiName=roi, DoseType='Max'),
                      'Average': dose.GetDoseStatistic(RoiName=roi,
                                                       DoseType='Average')}
        doses = dose.GetDoseAtRelativeVolumes(RoiName=roi,
                                              RelativeVolumes=volumes)
        for ii in range(len(volumes)):
            stats[roi][volume_names[ii]] = doses[ii]
    return stats


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
    """
    Hardcoded for now:
    pars = [
        1, PTV, MinDose, DoseLevel, 6200
        2, Rib, MaxDvh, DoseLevel, 3200
        3, SpinalCanal, MaxDvh, DoseLevel, 2080
        4, Heart, MaxDvh, DoseLevel, 2800
        5, Chestwall_L, MaxDvh, DoseLevel, 3000
        6, Lungs, MaxDvh, DoseLevel, 2000
    ]
    """
    start = time()

    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Set ROIs
    roi_names = ['PTV', 'Rib', 'SpinalCanal', 'Heart', 'Chestwall_L', 'Lungs']
    x0 =  [6200, 3200, 2080, 2800, 3000, 2000]
    dimensions = [(5500, 6200), (1600, 3200), (1040, 2080),
                  (1400, 2800), (1500, 3000), (1000, 2000)]
    
    # Get initial point
    print('Getting initial point')   
    y0 = objective(plan, beam_set, roi_names, x0)

    # Optimize
    print('\nStarting optimization')
    obj = lambda pars: objective(plan, beam_set, roi_names, pars)
    results = gp_minimize(obj, dimensions=dimensions, x0=x0, y0=y0,
                          random_state=0, n_calls=50, verbose=True)
    
    # need a callback to save progress if something goes wrong!
    
    # Save plan results
    print('\nSaving results')
    fpath = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\results\\'
    np.save(fpath + 'x.npy', results.x)
    np.save(fpath + 'x_iters.npy', results.x_iters)
    np.save(fpath + 'fun.npy', results.fun)
    np.save(fpath + 'func_vals.npy', results.func_vals)
    set_pars(plan, results.x)
    calc_plan(plan, beam_set);
    save_results(plan, roi_names, fpath)
    
    # Print total time
    total_time = time() - start
    print(f'Total time: {total_time} seconds')
 