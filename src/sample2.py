"""Sample treatment plans to discover correlations (patient 2)."""
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
        beam_set.NormalizeToPrescription(RoiName='PTV 4/7/20',
                                         DoseValue=6270.0, DoseVolume=99.0,
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
    volumes = [0.99, 0.98, 0.95, 0.9, 0.5, 0.1, 0.05, 0.02, 0.01]
    volume_names = ['D99', 'D98', 'D95', 'D90', 'D50', 'D10', 'D5', 'D2', 'D1']
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
    Objective terms = { Name : [{ idx, DoseLevel }] }
    Dose stats = { Name : { blah, blah, blah } }
    Clinical goals = { Name : [{ idx, AcceptanceLevel, GoalCriteria, GoalValue }] }
    Obj values = { Name : [{ FunctionType, DoseValue, PercentVolume, ResultValue }]}
    """    
    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Set ROIs
    roi_list = ['PTV 4/7/20', 'Lungs', 'SpinalCord (Thorax)',
                 'Esophagus', 'Heart']
    
    # Set objective function parameters
    par_dict = {
        'PTV 4/7/20': [{'idx': 1, 'DoseLevel': 7550, 'Range': [6271, 7550]}],
        'Lungs': [{'idx': 2, 'DoseLevel': 2000, 'Range': [0, 2000]}],
        'SpinalCord (Thorax)': [{'idx': 3, 'DoseLevel': 5000, 'Range': [0, 5000]}],
        'Esophagus': [{'idx': 4, 'DoseLevel': 6930, 'Range': [0, 6930]}],
        'Heart': [{'idx': 5, 'DoseLevel': 3500, 'Range': [0, 3500]}]
    }

    # Sample treatment plans
    n_success = 0
    results = []
    prev_pars = set()
    prev_pars.add(get_pars(par_dict))
    fpath = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\results\\'
    for ii in range(1000):
        print(f'Iteration: {ii}')
        if ii > 0:
            sample_pars(par_dict, prev_pars)
        print(get_pars(par_dict))
        set_pars(plan, par_dict)
        flag = calc_plan(plan, beam_set)
        n_success = n_success + 1 if flag == 0 else n_success
        print(f'Success: {flag}, n_success: {n_success}\n')
        if flag < 2:
            results.append([flag, copy.deepcopy(par_dict), get_stats(plan, roi_list),
                            get_goals(plan, roi_list), get_objs(plan)])
        else:
            results.append([flag, copy.deepcopy(par_dict)])
        np.save(fpath + 'samples_6_20.npy', results)
        if n_success >= 250:
            break
