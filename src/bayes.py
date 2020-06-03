"""Solve treatment plan with Bayesian optimization."""
from time import time

import numpy as np
from skopt import gp_minimize

import connect


# Eventually will need a better way to do objective, but for now maybe just
# need to stick with objective function terms indices...
# 0 SpinalCanal MaxDose
# 1 Esophagus MaxDose
# 2 Heart MaxDose
# 3 GreatVes MaxDose
# 4 Trachea MaxDose
# 5 Bronchus MaxDose
# 6 Rib MaxDose
# 7 Skin MaxDose
# 8 PTV MaxDose
# 9 PTV MinDose
# 10 Chestwall_L MaxDvh
# 11 Lungs MaxDvh
# 12 D2cm MaxDose
# 13 Pericardium MaxDose
# 14 R1 MaxDose
# 15 R2 MaxDose
# 16 D2cm_CD MaxDose
# 17 D2cm_CD2 MaxDose

# Need to make sure PTV MaxDose >= MinDose to avoid errors

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
    return score_plan(plan, roi_names) + [1e6, 0][norm]


def set_pars(plan, pars):
    """Set objective function parameters.

    Hardcoded for now:
    pars = [
        0 SpinalCanal, MaxDose, DoseLevel, 2080
        2 Heart MaxDose, DoseLevel, 2800
        6 Rib MaxDose, DoseLevel, 3200
        8 PTV MaxDose, DoseLevel, 6240
        9 PTV MinDose, DoseLevel, 4800
        10 Chestwall_L, MaxDvh, DoseLevel, 3000
        11 Lungs, MaxDvh, DoseLevel, 2000
        10 Chestwall_L, MaxDvh, PercentVolume, 1.5
        11 Lungs, MaxDvh, PercentVolume, 10
    ]

    Parameters
    ----------
    plan : connect.connect_cpytyon.PyScriptObject
        Current treatment plan.
    pars : list
        Objective function parameters.

    """
    funcs = plan.PlanOptimizations[0].Objective.ConstituentFunctions
    # roi_idx = [0, 2, 6, 8, 9, 10, 11]
    roi_idx = [0, 2, 6, 10, 11]
    for ii in range(len(roi_idx)):
        funcs[roi_idx[ii]].DoseFunctionParameters.DoseLevel = pars[ii]
    #funcs[10].DoseFunctionParameters.PercentVolume = pars[-2]
    #funcs[11].DoseFunctionParameters.PercentVolume = pars[-1]


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
    clinical_goals = plan.TreatmentCourse.EvaluationSetup.EvaluationFunctions
    for goal in clinical_goals:
        if goal.ForRegionOfInterest.Name in roi_names:
            criteria = goal.PlanningGoal.GoalCriteria
            level = goal.PlanningGoal.AcceptanceLevel
            try:
                value = goal.GetClinicalGoalValue()
                score += (-1)**(criteria == 'AtMost')*(level - value)/level
            except:
                return 0
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
    volumes = [0.99, 0.98, 0.95, 0.5, 0.05, 0.02, 0.01]
    volume_names = ['D99', 'D98', 'D95', 'D50', 'D5', 'D2', 'D1']
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
        0 SpinalCanal, MaxDose, DoseLevel, 2080
        2 Heart MaxDose, DoseLevel, 2800
        6 Rib MaxDose, DoseLevel, 3200
        8 PTV MaxDose, DoseLevel, 6240
        9 PTV MinDose, DoseLevel, 4800
        10 Chestwall_L, MaxDvh, DoseLevel, 3000
        11 Lungs, MaxDvh, DoseLevel, 2000
        10 Chestwall_L, MaxDvh, PercentVolume, 1.5
        11 Lungs, MaxDvh, PercentVolume, 10
    ]
    """
    start = time()

    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Set ROIs
    roi_names = ['SpinalCanal', 'Heart', 'Rib', 'PTV', 'Chestwall_L', 'Lungs']
    # x0 = [2080, 2800, 3200, 6240, 4800, 3000, 2000, 1.5, 10]
    x0 =  [2080, 2800, 3200, 3000, 2000] #, 1.5, 10]
    dimensions = [(1040, 2080), (1400, 2800), (1600, 3200), (1500, 3000),
                  (1000, 2000)]
                #  (0, 1.5), (0, 10)]
    
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
    calc_plan(plan, beam_set);
    save_results(plan, roi_names, fpath)
    
    # Print total time
    total_time = time() - start
    print(f'Total time: {total_time} seconds')
 