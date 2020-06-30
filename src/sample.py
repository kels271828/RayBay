"""Sample treatment plans on RayStation.

Constituent Functions
---------------------
The constituent function specifications are stored in a DataFrame with
columns corresponding to each function parameter: Roi, FunctionType,
DoseLevel, PercentVolume, EudParameterA, and Weight. The row index
within the DataFrame should correspond to the function index in the
RayStation objective function. Fixed parameters should be a single
value, tunable parameteres should be a list containing the minimum and
maximum values, and irrelevant parameters can be left blank. The
DataFrame can be built manually or from a CSV file using load_funcs().

Sampled Parameters
------------------
The sampled parameter values are stored in a DataFrame with columns
corresponding to each sample and function parameter: Sample, Term, Roi,
DoseLevel, PercentVolume, EudParameterA, and Weight. The term column
corresponds to the rows in the constituent function DataFrame.

Clinical Goals
--------------
The clinical goal specifications are stored in a DataFrame with columns
corresponding to each goal parameter: Roi, Type, GoalCriteria,
AcceptanceLevel, and ParameterValue. The DataFrame can be built
manually, based on the maximum parameter values in the funcs DataFrame
with init_goals(), or from a CSV file. The function get_results() is
currently able to evaluate MinDose, AverageDose, MaxDose, MinDvh, and
MaxDvh clinical goals. All other clinical goals are not evaluated, and
results are set to NaN.

Goal Results
------------
The clinical goal results are stored in a DataFrame with columns
Sample, Flag, and indices corresponding to the rows in the clinical
goals DataFrame.

Dose statistics
---------------
Dose statistics for a given list of regions of interest are stored in a
DataFrame with columms Sample, Roi, Min, Average, Max, D99, D98, D95,
D90, D50, D10, sD5, D2, and D1.

"""
import re

import numpy as np
import pandas as pd

import connect


def load_funcs(fpath):
    """Load constituent functions from CSV.

    Constituent function parameters should be specified by columns
    Roi, FunctionType, DoseLevel, PercentVolume, EudParameterA, and
    Weight. The row index within the CSV should correspond to the
    constituent function in the RayStation objective. Fixed parameters
    should be a single value, tunable parameteres should be a list
    containing the minimum and maximum values, and irrelevant
    parameters can be left blank.

    Parameters
    ----------
    fpath : str
        File path to CSV file.

    Returns
    -------
    pandas.DataFrame
        Constituent function specifications.

    """
    funcs = pd.read_csv(fpath).astype(object)
    for index, row in funcs.iterrows():
        for col in ['DoseLevel', 'PercentVolume', 'EudParameterA', 'Weight']:
            # Tunable parameters are read in as strings '[min, max]',
            # so we need to convert them back to a list of floats.
            if isinstance(row[col], str):
                pars = [float(par) for par
                        in re.findall(r'\d+\.\d+|\d+', row[col])]
                funcs.loc[index, col] = pars if len(pars) > 1 else pars[0]
    return funcs


def init_pars(funcs):
    """Initialize constituent function parameters.

    Constituent function parameters are specified by columns Sample,
    Term, Roi, DoseLevel, PercentVolume, EudParameterA, and Weight. The
    term columnn corresponds to the rows in the constituent function
    DataFrame.

    DoseLevel, PercentVolume, and EudParameterA are assigned to min or
    max value based on FunctionType. Weight is assigned to min value.

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
        'Term': index,
        'Roi': row['Roi'],
        'DoseLevel': get_par_bound(row['DoseLevel'], row['FunctionType']),
        'PercentVolume': get_par_bound(row['PercentVolume'],
                                       row['FunctionType']),
        'EudParameterA': get_par_bound(row['EudParameterA'],
                                       row['FunctionType']),
        'Weight': np.min(row['Weight'])
    } for index, row in funcs.iterrows()]
    columns = ['Sample', 'Term', 'Roi', 'DoseLevel', 'PercentVolume',
               'EudParameterA', 'Weight']
    return pd.DataFrame(data=data, columns=columns)


def get_par_bound(par, func_type):
    """Get min or max parameter value based on function type.

    Parameters
    ----------
    par : float or list
        Parameter value or boundaries.
    func_type : str
        Constituent function type.

    Returns
    -------
    float
        Parameter boundary value.

    """
    return np.max(par) if 'Max' in func_type else np.min(par)


def init_goals(funcs):
    """Initialize clinical goals based on constituent functions.

    Clinical goals are specified by columns Roi, Type, GoalCriteria,
    AcceptanceLevel, and ParameterValue. The function get_results() is
    currently able to evaluate MinDose, AverageDose, MaxDose, MinDvh,
    and MaxDvh clinical goals. All other clinical goals are not
    evaluated, and results are set to NaN.

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
        'Roi': row['Roi'],
        'Type': row['FunctionType'],
        'GoalCriteria': 'AtMost' if 'Max' in row['FunctionType'] else 'AtLeast',
        'AcceptanceLevel': get_par_bound(row['DoseLevel'],
                                         row['FunctionType']),
        'ParameterValue': get_par_bound(row['EudParameterA'],
                                        row['FunctionType'])
                          if 'Eud' in row['FunctionType'] else
                          get_par_bound(row['PercentVolume'],
                                        row['FunctionType'])
    } for index, row in funcs.iterrows()]
    columns = ['Roi', 'Type', 'GoalCriteria', 'AcceptanceLevel',
               'ParameterValue']
    return pd.DataFrame(data=data, columns=columns)


def init_results(goals):
    """Initialize clinical goal results.

    Clinical goal results are specified by Sample, Flag, and indices
    corresponding to the rows in the clinical goals DataFrame. Only
    MinDose, AverageDose, MaxDose, MinDvh, and MaxDvh are evaluated.
    All other clinical goals are set to NaN.

    Parameters
    ----------
    goals : pandas.DataFrame
        Clinical goal specifications.

    Returns
    -------
    pandas.DataFrame
        Clinical goal results.

    """
    columns = ['Sample', 'Flag'] + list(np.arange(len(goals)))
    return pd.DataFrame(columns=columns)


def init_stats():
    """Initialize dose statistics.

    Dose statistics are specified by Sample, Roi, Min, Average, Max,
    D99, D98, D95, D90, D50, D10, D5, D2, and D1.

    Returns
    -------
    pandas.DataFrame
        Dose statistics.

    """
    columns = ['Sample', 'Roi', 'Min', 'Average', 'Max', 'D99',
               'D98', 'D95', 'D90', 'D50', 'D10', 'D5', 'D2', 'D1']
    return pd.DataFrame(columns=columns)


def sample_pars(sample, funcs, pars):
    """Sample constituent function parameters.

    Constituent function parameters are specified by columns Sample,
    Term, Roi, DoseLevel, PercentVolume, EudParameterA, and Weight. The
    term columnn corresponds to the rows in the constituent function
    DataFrame.

    Parameters
    ----------
    sample : int
        Current sample number.
    funcs : pandas.DataFrame
        Constituent function specifications.
    pars : pandas.DataFrame
        Sampled constituent function parameters.

    Returns
    -------
    pandas.DataFrame
        Updated sampled constituent function parameters.

    """
    new_pars = []
    for index, row in funcs.iterrows():
        new_row = {'Sample': sample, 'Term': index, 'Roi': row['Roi']}
        new_row.update(sample_func_pars(row))
        new_pars.append(new_row)
    return pars.append(new_pars, ignore_index=True)


def sample_func_pars(func):
    """Sample constituent function parameters.

    Parameters
    ----------
    func : pandas.core.series.Series
        Row of constituent function specifications DataFrame.

    Returns
    -------
    pandas.core.series.Series
        Row of sampled constituent function parameter DataFrame.

    """
    pars = ['DoseLevel', 'PercentVolume', 'EudParameterA', 'Weight']
    new_row = {}
    for par in pars:
        if isinstance(func[par], list):
            low = func[par][0]
            high = func[par][1] + 1
            new_row[par] = np.random.randint(low, high)
        else:
            new_row[par] = func[par]
    return new_row


def set_pars(plan, pars):
    """Set objective function parameters.

    Constituent function parameters are specified by columns Sample,
    Term, Roi, DoseLevel, PercentVolume, EudParameterA, and Weight. The
    term columnn corresponds to the rows in the constituent function
    DataFrame.

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
        func = funcs[row['Term']].DoseFunctionParameters
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
        Region of interest for normalization.
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


def get_results(plan, sample, flag, goals, results):
    """Get clinical goal results.

    Clinical goals are specified by columns Roi, Type, GoalCriteria,
    AcceptanceLevel, and ParameterValue. Only MinDose, AverageDose,
    MaxDose, MinDvh, and MaxDvh are evaluated. All other clinical goals
    are set to NaN.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    sample : int
        Current sample number.
    flag : int
        0 = success, 1 = normalization failed, 2 = optimization failed
    goals : pandas.DataFrame
        Clinical goal specifications.
    results : pandas.DataFrame
        Clinical goal results.

    Returns
    -------
    pandas.DataFrame
        Updated clinical goal results.

    """
    dose = plan.TreatmentCourse.TotalDose
    new_results = {'Sample': sample, 'Flag': flag}
    if flag < 2:
        for index, row in goals.iterrows():
            new_results[index] = get_goal_result(dose, row)
    return results.append(new_results, ignore_index=True)


def get_goal_result(dose, goal):
    """Get result for given clinical goal.

    Only MinDose, AverageDose, MaxDose, MinDvh, and MaxDvh are
    evaluated. All other clinical goals are set to NaN.

    Parameters
    ----------
    dose : connect.connect_cpython.PyScriptObject
        Current treatment plan dose.
    goal : pandas.core.series.Series
        Row of clinical goal specification DataFrame.

    Returns
    -------
    float
        Clinical goal value.

    """
    if 'Dose' in goal['Type']:
        dose_type = re.findall('[A-Z][^A-Z]*', goal['Type'])[0]
        return dose.GetDoseStatistic(RoiName=goal['Roi'], DoseType=dose_type)
    elif 'Dvh' in goal['Type']:
        volume = 0.01*goal['ParameterValue']
        return dose.GetDoseAtRelativeVolumes(RoiName=goal['Roi'],
                                             RelativeVolumes=[volume])[0]
    else:
        return np.nan


def get_stats(plan, sample, roi_names, stats):
    """Get dose statistics for given regions of interest.

    Dose statistics are specified by Sample, Roi, Min, Average, Max,
    D99, D98, D95, D90, D50, D10, D5, D2, and D1.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    sample : int
        Current sample number.
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
    new_stats = []
    for roi in roi_names:
        new_row = {'Sample': sample, 'Roi': roi}
        new_row.update(get_roi_stats(dose, roi))
        new_stats.append(new_row)
    return stats.append(new_stats, ignore_index=True)


def get_roi_stats(dose, roi):
    """Get dose statistics for given region of interest.

    Dose statistics are specified by Sample, Roi, Min, Average, Max,
    D99, D98, D95, D90, D50, D10, D5, D2, and D1.

    Parameters
    ----------
    dose : conect.connect_cpython.PyScriptObject
        Current treatment plan dose.
    roi : str
        Region of interest to evaluate.

    Results
    -------
    dict
        Dose statistics.

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


def sample_plans(funcs, roi, dose, volume, goals=None, fpath='',
                 roi_names=None, max_iter=1000, n_success=100):
    """Sample treatment plans and save results.

    Results are saved after each iteration in case connection to
    RayStation times out.

    Parameters
    ----------
    funcs : pandas.DataFrame or str
        Constituent function specifications or path to CSV file.
    roi : str
        Region of interest for normalization.
    dose : float
        Dose value to normalize plans.
    volume : float
        Dose volume to normalize plans.
    goals : pandas.DataFrame or str, optional
        Clinical goal specifications or path to CSV file.
        If None, goals are based on constituent functions.
    fpath : str, optional
        Path to save results.
        If not specified, results are saved to the current directory.
    roi_names : iterable, optional
        Regions of interest to evaluate dose statistics.
        If None, based on regions of interest in clinical goals.
    max_iter : int, optional
        Maximum number of treatment plans to sample.
    n_success : int, optional
        Number of successfull treatment plans to sample.

    Returns
    -------
    pandas.DataFrame
        Sampled constituent function parameters.
    pandas.DataFrame
        Clinical goal results.
    pandas.DataFrame
        Dose statistics.

    """
    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Define functions and goals
    funcs, goals, pars, results, stats = init_prob(funcs, goals, roi_names)

    # Sample treatment plans
    count = 0
    for ii in range(max_iter):
        print(f'Iteration: {ii}', end='')
        if ii > 0:
            pars = sample_pars(ii, funcs, pars)
            pars.to_pickle(fpath + 'pars.npy')
        set_pars(plan, pars)
        flag = calc_plan(plan, beam_set, roi, dose, volume)
        count = count + 1 if flag == 0 else count
        print(f', Flag: {flag}, Successes: {count}')
        results = get_results(plan, ii, flag, goals, results)
        results.to_pickle(fpath + 'results.npy')
        if flag < 2:
            stats = get_stats(plan, ii, roi_names, stats)
            stats.to_pickle(fpath + 'stats.npy')
        if count == n_success:
            break
    return pars, results, stats


def init_prob(funcs, goals, roi_names):
    """Initialize treatment plan sampling structures.

    Parameters
    ----------
    funcs : pandas.DataFrame or str
        Constituent function specifications or path to CSV file.
    goals : pandas.DataFrame or str, optional
        Clinical goal specifications or path to CSV file.
    roi_names : iterable, optional
        Regions of interest to evaluate dose statistics.

    Returns
    -------
    pandas.DataFrame
        Constituent function specifications.
    pandas.DataFrame
        Clinical goal specifications.
    pandas.DataFrame
        Sampled constituent function parameters.
    pandas.DataFrame
        Clinical goal results.
    pandas.DataFrame
        Dose statistics.

    """
    if isinstance(funcs, str):
        funcs = load_funcs(funcs)
    if goals is None:
        goals = init_goals(funcs)
    elif isinstance(goals, str):
        goals = pd.read_csv(goals)
    if roi_names is None:
        roi_names = set(goals['Roi'])
    pars = init_pars(funcs)
    results = init_results(goals)
    stats = init_stats()
    return funcs, goal, pars, results, stats


def grid_search(funcs, roi, dose, volume, goals=None, fpath='',
                 roi_names=None):
    """Perform grid search over treatment plans and save results.

    Results are saved after each iteration in case connection to
    RayStation times out.

    Parameters
    ----------
    funcs : pandas.DataFrame or str
        Constituent function specifications or path to CSV file.
    roi : str
        Region of interest for normalization.
    dose : float
        Dose value to normalize plans.
    volume : float
        Dose volume to normalize plans.
    goals : pandas.DataFrame or str, optional
        Clinical goal specifications or path to CSV file.
        If None, goals are based on constituent functions.
    fpath : str, optional
        Path to save results.
        If not specified, results are saved to the current directory.
    roi_names : iterable, optional
        Regions of interest to evaluate dose statistics.
        If None, based on regions of interest in clinical goals.

    Returns
    -------
    pandas.DataFrame
        Sampled constituent function parameters.
    pandas.DataFrame
        Clinical goal results.
    pandas.DataFrame
        Dose statistics.

    """
    # Get RayStation Objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Define functions and goals
    funcs, goals, pars, results, stats = init_prob(funcs, goals, roi_names)

    # Sample treatment plans
    # Divide into how many grid points???
    # Log scale or linear scale???
    for ii in range(max_iter):
        print(f'Iteration: {ii}', end='')
        if ii > 0:
            pars = sample_pars(ii, funcs, pars)
            pars.to_pickle(fpath + 'pars.npy')
        set_pars(plan, pars)
        flag = calc_plan(plan, beam_set, roi, dose, volume)
        count = count + 1 if flag == 0 else count
        print(f', Flag: {flag}, Successes: {count}')
        results = get_results(plan, ii, flag, goals, results)
        results.to_pickle(fpath + 'results.npy')
        if flag < 2:
            stats = get_stats(plan, ii, roi_names, stats)
            stats.to_pickle(fpath + 'stats.npy')
    return pars, results, stats









