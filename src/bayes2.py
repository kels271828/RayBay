"""Solve treatment plan with Bayesian optimization (patient 2)."""
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
        return 1e6


def set_pars(plan, pars):
    """Set objective function parameters.

    Hardcoded for now:
    pars = [
        0, PTV, MinDose, DoseLevel, 6270
        1, PTV, MaxDose, DoseLevel, 7550
        2, Lungs, MaxDvh, DoseLevel, 2000
        3, SpinalCord, MaxDose, DoseLevel, 5000
        4, Esophagus, MaxDose, DoseLevel, 6930
        5, Heart, MaxEud, DoseLevel, 3500
    ]

    Parameters
    ----------
    plan : connect.connect_cpytyon.PyScriptObject
        Current treatment plan.
    pars : list
        Objective function parameters.

    """
    funcs = plan.PlanOptimizations[0].Objective.ConstituentFunctions
    roi_idx = [1, 2, 3, 4, 5]
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
        print('Could not optimize')
        return False

    # Normalize plan
    try:
        beam_set.NormalizeToPrescription(RoiName='PTV 4/7/20',
                                         DoseValue=6270.0, DoseVolume=99.0,
                                         PrescriptionType='DoseAtVolume')
        return True
    except:
        print('Could not normalize')
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

    """
    score = 0
    vals = get_objs(plan)
    for roi in roi_names:
        level = vals[roi][-1]['DoseValue']
        value = vals[roi][-1]['ResultValue']
        if roi == 'PTV 4/7/20':
            if value < level:
                score += 1/(level - value)
            else:
                return 1e6
        else:
            score += 100*(value - level)/level
    return score


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
    np.save(fpath + 'objs.npy', get_objs(plan))


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
    values['PTV 4/7/20'] = [{
        'FunctionType': 'MinDose',
        'DoseValue': 6270,
        'PercentVolume': 0,
        'ResultValue': dose.GetDoseStatistic(RoiName='PTV 4/7/20', DoseType='Min')
    }, {
        'FunctionType': 'MaxDose',
        'DoseValue': 7550,
        'PercentVolume': 0,
        'ResultValue': dose.GetDoseStatistic(RoiName='PTV 4/7/20', DoseType='Max')
    }]
    values['Lungs'] = [{
        'FunctionType': 'MaxDvh',
        'DoseValue': 2000,
        'PercentVolume': 35,
        'ResultValue': dose.GetDoseAtRelativeVolumes(RoiName='Lungs',
            RelativeVolumes=[0.35])[0]
    }]
    values['SpinalCord (Thorax)'] = [{
        'FunctionType': 'MaxDose',
        'DoseValue': 5000,
        'PercentVolume': 0,
        'ResultValue': dose.GetDoseStatistic(RoiName='SpinalCord (Thorax)', DoseType='Max')
    }]
    values['Esophagus'] = [{
        'FunctionType': 'MaxDose',
        'DoseValue': 6930,
        'PercentVolume': 0,
        'ResultValue': dose.GetDoseStatistic(RoiName='Esophagus', DoseType='Max')
    }]
    values['Heart'] = [{
        'FunctionType': 'MaxAvg',
        'DoseValue': 3500,
        'PercentVolume': 0,
        'ResultValue': dose.GetDoseStatistic(RoiName='Heart', DoseType='Average')
    }]
    return values


if __name__ == '__main__':
    """
    Hardcoded for now:
    pars = [
        0, PTV, MinDose, DoseLevel, 6270
        1, PTV, MaxDose, DoseLevel, 7550
        2, Lungs, MaxDvh, DoseLevel, 2000
        3, SpinalCord, MaxDose, DoseLevel, 5000
        4, Esophagus, MaxDose, DoseLevel, 6930
        5, Heart, MaxEud, DoseLevel, 3500
    ]
    """
    start = time()

    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Set ROIs
    roi_names = ['PTV 4/7/20', 'Lungs', 'SpinalCord (Thorax)',
                 'Esophagus', 'Heart']
    x0 =  [6271, 2000, 5000, 6930, 3500]
    dimensions = [(6271, 6590), (1875, 3125), (3750, 6250),
                  (5197.5, 8662.5), (2625, 4375)]
    #dimensions = [(6910, 7550), (1000, 2000), (2500, 5000),
    #              (3465, 6930), (1750, 3500)]
    #dimensions = [(6590., 7550.), (500., 2000.), (1250., 5000.),
    #              (1732.5, 6930.), (875., 3500.)]
    
    # Get initial point
    print('Getting initial point')   
    y0 = objective(plan, beam_set, roi_names, x0)
    #fpath = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\results\\'
    #temp = np.load(fpath + '6_16_quarter\\x_iters.npy')
    #x0 = [list(temp[ii]) for ii in range(temp.shape[0])]
    #y0 = list(np.load(fpath + '6_16_quarter\\func_vals.npy'))

    # Optimize
    print('\nStarting optimization')
    obj = lambda pars: objective(plan, beam_set, roi_names, pars)
    results = gp_minimize(obj, dimensions=dimensions, x0=x0, y0=y0,
                          random_state=0, n_calls=50, verbose=True)
    #results = gp_minimize(obj, dimensions=dimensions, x0=x0, y0=y0,
    #                      random_state=0, n_calls=50, verbose=True,
    #                      n_random_starts=0)
    
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
 