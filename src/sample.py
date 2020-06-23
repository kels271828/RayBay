"""Sample treatment plans to discover correlations.

The constituent functions DataFrame should have columns 'Roi',
'DoseLevel', 'PercentVolume', 'EudParameterA', and 'Weight'. Rows
should be in the order that they appear in the objective function.
Fixed parameters should be a single value, and tunable parameters
should be a list containting the minimum and maximum allowable values.

The parameter DataFrame has columns 'Sample', 'Term', 'Roi',
'DoseLevel', 'PercentVolume', 'EudParameterA', and 'Weight'.

The clinical goals DataFrame has columns 'Roi', 'Type', 'GoalCriteria',
'AcceptanceLevel', and 'ParameterValue'.

The results DataFrame has columns 'Sample', 'Flag', and indices
corresponding to each clinical goal.

The stats DataFrame has columns 'Sample', 'Roi', 'Min', 'Average',
'Max', 'D99', 'D98', 'D95', 'D90', 'D50', 'D10', 'D5', 'D2, and 'D1'.

"""
import copy

import numpy as np
import pandas as pd

import connect


def init_funcs(fpath):
    """Initialize constituent functions from CSV.
    
    Parameters
    ----------
    fpath : str
        File path to CSV file.
    
    Returns
    -------
    pandas.DataFrame
        Constituent function specifications.
    
    """
    funcs = pd.read_csv(fpath)
    funcs = funcs.astype(object)
    col_names = ['DoseLevel', 'PercentVolume', 'EudParameterA', 'Weight']
    for index, row in funcs.iterrows():
        for col in col_names:
            if isinstance(row[col], str):
                temp = [float(s) for s in re.findall(r'\d+\.\d+|\d+', row[col])]
                funcs.loc[index, col] = (temp, temp[0])[len(temp) == 1]
    return funcs


def init_pars(funcs):
    """Initialize constituent function parameters to maximum values.
    
    Parameters
    ----------
    funcs : pandas.DataFrame
        Constituent function specifications.
    
    Returns
    -------
    pandas.DataFrame
        Sampled constituent function parameters.
    
    """
    data = [{
        'Sample': 0,
        'Term': ii,
        'Roi': funcs.iloc[ii]['Roi'],
        'DoseLevel': np.max(funcs.iloc[ii]['DoseLevel']),
        'PercentVolume': np.max(funcs.iloc[ii]['PercentVolume']),
        'EudParameterA': np.max(funcs.iloc[ii]['EudParameterA']),
        'Weight': np.max(funcs.iloc[ii]['Weight'])
    } for ii in range(6)]
    columns = ['Sample', 'Term', 'Roi', 'DoseLevel', 'PercentVolume',
               'EudParameterA', 'Weight']
    return pd.DataFrame(data=data, columns=columns)


def init_goals(funcs):
    """Initialize clinical goals based on constituent functions.
    
    Parameters
    ----------
    funcs : pandas.DataFrame
        Constituent function specifications.
    
    Returns
    -------
    pandas.DataFrame
        Clinical goal specifications.
    
    """
    data = [{
        'Roi': funcs.iloc[ii]['Roi'],
        'Type': funcs.iloc[ii]['FunctionType'],
        'GoalCriteria': ('AtLeast', 'AtMost')\
                        ['Max' in funcs.iloc[ii]['FunctionType']],
        'AcceptanceLevel': np.max(funcs.iloc[ii]['DoseLevel']),
        'ParameterValue': (funcs.iloc[ii]['PercentVolume'],
                           funcs.iloc[ii]['EudParameterA'])\
                          ['Eud' in funcs.iloc[ii]['FunctionType']]
    } for ii in range(6)]
    columns = ['Roi', 'Type', 'GoalCriteria', 'AcceptanceLevel',
               'ParameterValue']
    return pd.DataFrame(data=data, columns=columns)


def init_results(goals):
    """Initialize clinical goal results.

    Parameters
    ----------
    goals : pandas.DataFrame
        Clinical goal specifications.

    Returns
    -------
    pandas.DataFrame
        Clinical goal results.

    """
    columns = ['Sample', 'Flag']
    columns += list(np.arange(len(goals)))
    return pd.DataFrame(columns=columns)


def init_stats():
    """Initialize dose statistics.

    Returns
    -------
    pandas.DataFrame
        Dose statistics.

    """
    columns = ['Sample', 'Roi', 'Min', 'Average', 'Max', 'D99',
               'D98', 'D95', 'D90', 'D50', 'D10', 'D5', 'D2', 'D1']
    return pd.DataFrame(columns=columns)


def sample_pars(funcs, pars):
    """Sample constituent function parameters.
    
    Parameters
    ----------
    funcs : pandas.DataFrame
        Constituent function specifications.
    pars : pandas.DataFrame
        Sampled constituent function parameters.
    
    Returns
    -------
    pandas.DataFrame
        Updated sampled constituent function parameters.
    
    """
    sample = _get_sample_num(pars)
    new_pars = []
    for idx, row in funcs.iterrows():
        new_row = {'Sample': sample, 'Term': idx, 'Roi': row['Roi']}
        new_row.update(_sample_func_pars, row)
        new_pars.append(new_row)
    return pars.append(new_pars, ignore_index=True)


def _sample_func_pars(func):
    """Sample constituent function parameters."""
    pars = ['DoseLevel', 'PercentVolume', 'EudParameterA', 'Weight']
    for par in pars:
        if isinstance(func[par], list):
            low = func[par][0]
            high = func[par][1] + 1
            new_row[par] = np.random.randint(low, high)
        else:
            new_row[par] = func[par]
    return new_row


def _get_sample_num(df):
    """Get current sample number."""
    if len(df) == 0:
        return 0
    else:
        return max(df['Sample']) + 1


def set_pars(plan, pars):
    """Set objective function parameters.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.        
    pars : pandas.DataFrame
        Constituent function parameters.

    Returns
    -------
    None.

    """
    funcs = plan.PlanOptimizations[0].Objective.ConstituentFunctions
    for _, row in pars[pars['Sample'] == max(pars['Sample'])].iterrows():
        func = funcs[row['Term']].DoseFunctionParameteres
        func.DoseLevel = row['DoseLevel']
        func.Weight = row['Weight']
        if func.FunctionType == 'MaxEud':
            func.EudParameterA = row['EudParameterA']
        else:
            func.PercentVolume = row['PercentVolume']


def calc_plan(plan, beam_set, roi, dose, volume):
    """Calculate and normalize treatment plan.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    beam_set : connect.connect_cpython.PyScriptObject
        Current beam set.
    roi : str
        ROI to normalize plan.
    dose : float
        Dose value to normalize plan.
    volume : float
        Dose volume to normalize plan.
        
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
        beam_set.NormalizeToPrescription(RoiName=roi, DoseValue=dose, 
                                         DoseVolume=volume,
                                         PrescriptionType='DoseAtVolume')
        return 0
    except:
        return 1


def get_results(plan, flag, goals, results):
    """Get clinical goal results.

    Parameters
    ----------
    flag : int
        0 = success, 1 = normalization failed, 2 = optimization failed
    pandas.DataFrame
        Clinical goal specifications.

    Returns
    -------
    pandas.DataFrame
        Clinical goal results.

    """
    dose = plan.TreatmentCourse.TotalDose
    sample = _get_sample_num(results)
    new_results = {'Sample': sample, 'Flag': flag}
    for index, row in goals.iterrows():
        new_results[index] = _get_roi_result(dose, row)
    return results.append(new_results, ignore_index=True)


def _get_roi_result(dose, goal):
    """Get result for given clinical goal."""
    if 'Dose' in goal['Type']:
        dose = re.findall('[A-Z][^A-Z]*', goal['Type'])[0]
        return dose.GetDoseStatistic(RoiName=goal['Roi'], DoseType=dose)
    elif 'Dvh' in goal['Type']:
        volume = 0.01*row['ParameterValue']
        return dose.GetDoseAtRelativeVolumes(RoiName=goal['Roi'],
                                             RelativeVolumes=[volume])[0]
    else:
        return np.nan


def get_stats(plan, roi_names, stats):
    """Get dose statistics for given regions of interest.
    
    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    roi_names : iterable
        Regions of interest to evaluate.
    stats : pandas.DataFrame
        Dose statistics.

    Returns
    -------
    pandas.DataFrame
        Updated dose statistics.

    """
    dose = plan.TreatmentCourse.TotalDose
    sample = _get_sample_num(stats)
    new_stats = []
    for roi in roi_list:
        new_row = {'Sample': sample, 'Roi': roi}
        new_row.update(_get_roi_stats(dose, roi))
        new_stats.append(new_row)
    return stats.append(new_stats, ignore_index=True)


def _get_roi_stats(dose, roi):
    """Get dose statistics for given region of interest."""
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
    roi_list = ['SpinalCanal', 'Heart', 'Rib', 'PTV', 'Chestwall_L', 'Lungs']
    
    # Set objective function parameters
    par_dict = {
        'PTV': [{'idx': 1, 'DoseLevel': 6200, 'Range': [4801, 6200]}],
        'Rib': [{'idx': 2, 'DoseLevel': 3200, 'Range': [0, 3200]}],
        'SpinalCanal': [{'idx': 3, 'DoseLevel': 2080, 'Range': [0, 2080]}],
        'Heart': [{'idx': 4, 'DoseLevel': 2800, 'Range': [0, 2800]}],
        'Chestwall_L': [{'idx': 5, 'DoseLevel': 3000, 'Range': [0, 3000]}],
        'Lungs': [{'idx': 6, 'DoseLevel': 2000, 'Range': [0, 2000]}]
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
        np.save(fpath + 'samples_6_11.npy', results)
        if n_success >= 250:
            break
