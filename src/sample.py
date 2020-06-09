"""Sample treatment plans to discover correlations."""
import copy

import numpy as np

import connect


def sample_pars(par_dict, prev_pars):
    """Sample objective function parameters with no repeats.

    Parameters should have the form

    { Name : [{ idx, DoseLevel }] }

    where idx is the index of the objective function term.

    Parameters
    ----------
    par_dict : dict
        Objective function parameters.
    prev_pars : set
        Results from previous samples.

    Returns
    -------
    None.

    """
    while get_pars(par_dict) in prev_pars:
        for roi in par_dict:
            for term in par_dict[roi]:
                term['DoseLevel'] = np.random.randint(term['Range'][0],
                                                      term['Range'][1])
    prev_pars.add(get_pars(par_dict))


def get_pars(par_dict):
    """Get tuple of parameters from dictionary.

    Parameters
    ----------
    par_dict : dict
        Objective function parameters.

    Returns
    -------
    par_tuple : tuple
        Objective function parameters.

    """
    return tuple(term['DoseLevel'] for roi in par_dict for term in par_dict[roi])


def set_pars(plan, par_dict):
    """Set objective function parameters.
    
    Parameters should have the form
    
    { Name : [{ idx, DoseLevel }] }

    where idx is the index of the objective function term.   

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.        
    par_dict : dict
        Objective function parameters.

    Returns
    -------
    None.

    """
    funcs = plan.PlanOptimizations[0].Objective.ConstituentFunctions
    for roi in par_dict:
        for term in par_dict[roi]:
            funcs[term['idx']].DoseFunctionParameters.DoseLevel = term['DoseLevel']


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
    int
        0 = success, 1 = normalization failed, 2 = optimization failed  

    """
    # Calculate plan
    plan.PlanOptimizations[0].ResetOptimization()
    try:
        plan.PlanOptimizations[0].RunOptimization()
    except:
        return 2

    # Normalize plan
    try:
        beam_set.NormalizeToPrescription(RoiName='PTV', DoseValue=4800.0,
                                         DoseVolume=95.0,
                                         PrescriptionType='DoseAtVolume')
        return 0
    except:
        return 1


def get_stats(plan, roi_list):
    """Get dose statistics for given regions of interest.
    
    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    roi_list : list
        Regions of interest to evaluate.
    Returns
    -------
    dict
        Dose statistics for given regions of interest.
    """
    stats = {}
    dose = plan.TreatmentCourse.TotalDose
    for roi in roi_list:
        stats[roi] = get_roi_stats(dose, roi)
    return stats


def get_roi_stats(dose, roi):
    """Get dose statistics for given region of interest.

    Parameters
    ----------
    dose : connect.connect_cpython.PyScriptObject
        Total dose for current treatment plan.
    roi : str
        Region of interest to evaluate.

    Returns
    -------
    dict
        Dose statistics for given region of interest.

    """
    volumes = [0.99, 0.98, 0.95, 0.5, 0.05, 0.02, 0.01]
    volume_names = ['D99', 'D98', 'D95', 'D50', 'D5', 'D2', 'D1']
    stats = {'Min': dose.GetDoseStatistic(RoiName=roi, DoseType='Min'),
             'Max': dose.GetDoseStatistic(RoiName=roi, DoseType='Max'),
             'Average': dose.GetDoseStatistic(RoiName=roi, DoseType='Average')}
    doses = dose.GetDoseAtRelativeVolumes(RoiName=roi,
                                          RelativeVolumes=volumes)
    for ii in range(len(volumes)):
        stats[volume_names[ii]] = doses[ii]
    return stats


def get_goals(plan, roi_list):
    """Get clinical goals for given regions of interest.
    
    Results have the form
    
    { Name : [{ idx, AcceptanceLevel, GoalCriteria, GoalValue }] }
    
    where idx is the index of the evaluation function term.
    
    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    roi_list : list
        Regions of interest to evaluate.
        
    Returns
    -------
    dict
        Clinical goals for given regions of interest.
    
    """
    # Initialize results
    results = {}
    for roi in roi_list:
        results[roi] = []
        
    # Get clinical goals    
    goals = plan.TreatmentCourse.EvaluationSetup.EvaluationFunctions
    for ii in range(len(goals)):
        roi = goals[ii].ForRegionOfInterest.Name
        if roi in roi_list:
            level = goals[ii].PlanningGoal.AcceptanceLevel
            criteria = goals[ii].PlanningGoal.GoalCriteria
            try:
                value = goals[ii].GetClinicalGoalValue()
            except:
                value = -1
            results[roi].append({'AcceptanceLevel': level, 'GoalValue': value,
                                 'GoalCriteria': criteria, 'idx': ii})
    return results


if __name__ == '__main__':
    """
    Objective function terms:
     *  0, SpinalCanal, MaxDose, 2080
        1, Esophagus, MaxDose
     *  2, Heart, MaxDose, 2800
        3, GreatVes, MaxDose
        4, Trachea, MaxDose
        5, Bronchus, MaxDose
     *  6, Rib, MaxDose, 3200 
        7, Skin, MaxDose
     *  8, PTV, MaxDose, 6240
     *  9, PTV, MinDose, 4800
     * 10, Chestwall_L, MaxDvh, 3000, 1.5
     * 11, Lungs, MaxDvh, 2000, 10
       12, D2cm, MaxDose
       13, Pericardium, MaxDose
       14, R1, MaxDose
       15, R2, MaxDose
       16, D2cm_CD, MaxDose
       17, D2cm_CD2, MaxDose

    Objective terms = { Name : [{ idx, DoseLevel }] }
    Dose stats = { Name : { blah, blah, blah } }
    Clinical goals = { Name : [{ idx, AcceptanceLevel, GoalCriteria, GoalValue }] }
    """    
    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Set ROIs
    roi_list = ['SpinalCanal', 'Heart', 'Rib', 'PTV', 'Chestwall_L', 'Lungs']
    
    # Set objective function parameters
    par_dict = {
        'SpinalCanal': [{'idx': 0, 'DoseLevel': 2080, 'Range': [0, 2080]}],
        'Heart': [{'idx': 2, 'DoseLevel': 2800, 'Range': [0, 2800]}],
        'Rib': [{'idx': 6, 'DoseLevel': 3200, 'Range': [0, 3200]}],
        #'PTV': [{'idx': 8, 'DoseLevel': 6240, 'Range': [4801, 6240]}],
        'Chestwall_L': [{'idx': 10, 'DoseLevel': 3000, 'Range': [0, 3000]}],
        'Lungs': [{'idx': 11, 'DoseLevel': 2000, 'Range': [0, 2000]}]
    }

    # Sample treatment plans
    results = []
    prev_pars = set()
    prev_pars.add(get_pars(par_dict))
    fpath = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\results\\'
    for ii in range(100):
        print(f'Iteration: {ii}')
        if ii > 0:
            sample_pars(par_dict, prev_pars)
        print(get_pars(par_dict))
        set_pars(plan, par_dict)
        flag = calc_plan(plan, beam_set)
        print(f'Success: {flag}\n')
        if flag == 0:
            results.append([copy.deepcopy(par_dict), get_stats(plan, roi_list),
                            get_goals(plan, roi_list)])
        else:
            results.append([copy.deepcopy(par_dict), flag])
        np.save(fpath + 'samples_6_8.npy', results)
 