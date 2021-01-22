"""RayStation treatment planning with Bayesian optimization.

TODO:
* Add 1D grid_search function

"""
import pickle
import re
from time import time

import numpy as np
import skopt

import connect
import raybay


def get_plan(funcs, norm, goals=None, solver='gp_minimize', n_calls=25,
             random_state=None, n_initial_points=10, verbose=True):
    """Hyperparameter optimization for RayStation treatment planning.

    Hyperparameter optimization for RayStation treatment planning using
    the following solvers from scikit-optimize:

        - `gp_minimize`: Bayesian optimization using Gaussian processes.
        - `forest_minimize`: Sequential optimization using decision
           trees.
        - `dummy_minimize`: Random search by uniform sampling within the
           given bounds.

    For more details about scikit-optimize, refer to
    https://scikit-optimize.github.io/stable/index.html

    Parameters
    ----------
    funcs : str
        Path to CSV with constituent function specifications.
    norm : (str, float, float)
        Region of interest, dose, and volume used for normalization.
    goals : pandas.DataFrame or str, optional
        Path to CSV with clinical goal specifications.
        If None, goals are assigned based on constituent functions.
    solver : {'gp_minimize', 'forest_minimize', 'dummy_minimize'}, optional
        Name of scikit-optimize solver to use.
    n_calls : int, optional
        Number of calls to objective.
    random_state : int, optional
        Set random state for reproducible results.
    n_initial_points : int, optional
        Number of random function evaluations before function
        approximation.
    verbose : bool, optional
        Control the verbosity of the solver.

    Returns
    -------
    raybay.RaybayResult
        RayStation treatment plan results.

    """
    # Get RayStation objects
    patient = connect.get_current('Patient')
    case = connect.get_current('Case')
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Initialize result object
    result = raybay.RaybayResult(patient.Name, case.CaseName, plan.Name, funcs,
                                 norm, goals, solver)

    # Optimize
    def obj(pars):
        return objective(plan, beam_set, result, funcs[:-9], pars)
    checkpoint_path = funcs[:-9] + 'checkpoint.pkl'
    checkpoint_saver = skopt.callbacks.CheckpointSaver(checkpoint_path,
                                                       store_objective=False)
    start_time = time()
    if solver == 'forest_minimize':
        result.opt_result = skopt.forest_minimize(
            obj,
            dimensions=get_dims(result.func_df),
            n_calls=n_calls,
            n_initial_points=n_initial_points,
            random_state=random_state,
            verbose=verbose,
            callback=[checkpoint_saver])
    elif solver == 'dummy_minimize':
        result.opt_result = skopt.dummy_minimize(
            obj,
            dimensions=get_dims(result.func_df),
            n_calls=n_calls,
            random_state=random_state,
            verbose=verbose,
            callback=[checkpoint_saver])
    else:
        result.opt_result = skopt.gp_minimize(
            obj,
            dimensions=get_dims(result.func_df),
            n_calls=n_calls,
            n_initial_points=n_initial_points,
            random_state=random_state,
            verbose=verbose,
            callback=[checkpoint_saver])
    result.opt_result.specs['args']['func'] = 'local'  # remove local func
    result.time = time() - start_time                  # to allow pickling

    # Get optimal dose-volume histogram
    set_pars(plan, result.func_df, result.opt_result.x)
    calc_plan(plan, beam_set, result.norm)
    result.dvh_dict = get_dvh(result.roi_list)

    return result


def objective(plan, beam_set, result, repo_path, pars):
    """Objective function for hyperparameter optimization.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    beam_set : connect.connect_cpython.PyScriptObject
        Current beam set.
    result : raybay.RaybayResult
        RayStation treatment plan results.
    repo_path : str
        Path to save goal checkpoints.
    pars : list
        Constituent function parameters.

    Returns
    -------
    float
        Treatment plan score.

    """
    set_pars(plan, result.func_df, pars)
    flag = calc_plan(plan, beam_set, result.norm)
    result.flag_list.append(flag)
    print(f'Flag: {flag}', flush=True)
    return get_score(plan, result.goal_df, result.norm, flag, result.goal_dict,
                     repo_path)


def set_pars(plan, func_df, pars):
    """Set objective function parameters.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    func_df : pandas.DataFrame
        Constituent function specifications.
    pars : list
        Constituent function parameters.

    Returns
    -------
    None.

    """
    count = 0
    const_funcs = plan.PlanOptimizations[0].Objective.ConstituentFunctions
    for index, row in func_df.iterrows():
        func = const_funcs[index].DoseFunctionParameters
        if isinstance(row['DoseLevel'], list):
            func.DoseLevel = pars[count]
            count += 1
        else:
            func.DoseLevel = row['DoseLevel']
        if 'Eud' in func.FunctionType:
            func.EudParameterA = row['EudParameterA']
        elif isinstance(row['PercentVolume'], list):
            func.PercentVolume = pars[count]
            count += 1
        else:
            func.PercentVolume = row['PercentVolume']
        if isinstance(row['Weight'], list):
            func.Weight = pars[count]
            count += 1
        else:
            func.Weight = row['Weight']


def calc_plan(plan, beam_set, norm):
    """Calculate and normalize treatment plan.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    beam_set : connect.connect_cpython.PyScriptObject
        Current beam set.
    norm : (str, float, float)
        Region of interest, dose, and volume used for normalization.

    Returns
    -------
    int
        RayStation exit status:
        - 0: success
        - 1: normalization failed
        - 2: optimization failed

    """
    # Calculate plan
    plan.PlanOptimizations[0].ResetOptimization()
    try:
        plan.PlanOptimizations[0].RunOptimization()
    except:
        return 2

    # Normalize plan
    try:
        beam_set.NormalizeToPrescription(
            RoiName=norm[0],
            DoseValue=norm[1],
            DoseVolume=norm[2],
            PrescriptionType='DoseAtVolume')
        return 0
    except:
        return 1


def get_score(plan, goal_df, norm, flag, goal_dict, repo_path):
    """Calculate treatment plan score.

    The treatment plan score is a linear combination of the relative
    difference between goal values and goal results. Returns 1e6 if
    RayStation optimization failed.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    goal_df : pandas.DataFrame
        Clinical goal specifications.
    norm : (str, float, float)
        Region of interest, dose, and volume used for normalization.
    flag : int
        RayStation exit status.
    goal_dict : dict
        Clinical goal results.
    repo_path : str
        Path to save goal checkpoints

    Returns
    -------
    float
        Treatment plan score.

    """
    if flag == 2:
        return 1e6
    results = get_results(plan, goal_df)
    scale = get_scale(goal_df, norm, results) if flag == 1 else 1.0
    score = 0
    for index, row in goal_df.iterrows():
        value = scale*results[index]
        goal_dict[index].append(value)
        score += -row['Weight']*raybay.get_term(value, row['AcceptanceLevel'],
                                                row['Type'], row['Shape'])
    with open(repo_path + 'goal_dict.pkl', 'wb') as fp:
        pickle.dump(goal_dict, fp)
    return score


def get_results(plan, goal_df):
    """Get clinical goal results.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    goal_df : pandas.DataFrame
        Clinical goal specifications.

    Returns
    -------
    dict
        Clinical goal results.

    """
    dose = plan.TreatmentCourse.TotalDose
    results = {}
    for index, row in goal_df.iterrows():
        results[index] = get_value(dose, row)
    return results


def get_value(dose, goal):
    """Get clinical goal value.

    Currently able to evaluate MinDose, AverageDose, MaxDose, MinDvh,
    and MaxDvh clinical goal values. All other clinical goal values are
    returned as NaN.

    Parameters
    ----------
    dose : connect.connect_cpython.PyScriptObject
        Current treatment plan dose.
    goal : pandas.core.series.series
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


def get_scale(goal_df, norm, results):
    """Get normalization scale factor.

    Parameters
    ----------
    goal_df : pandas.DataFrame
        Clinical goal specifications.
    norm : (str, float, float)
        Region of interest, dose, and volume used for normalization.
    results : dict
        Clinical goal results for current iteration.

    Returns
    -------
    float
        Normalization scale factor.

    """
    index = goal_df.index[(goal_df['Roi'] == norm[0]) &
                          (goal_df['AcceptanceLevel'] == norm[1]) &
                          (goal_df['ParameterValue'] == norm[2])].tolist()[0]
    return norm[1]/results[index]


def get_dims(func_df):
    """Get constituent function parameter dimensions.

    Parameters
    ----------
    func_df : pandas.DataFrame
        Constituent function specifications.

    Returns
    -------
    list
        Constituent function parameter dimensions.

    """
    dimensions = []
    for _, row in func_df.iterrows():
        for par in ['DoseLevel', 'PercentVolume', 'Weight']:
            if isinstance(row[par], list):
                dimensions.append(row[par])
    return dimensions


def get_dvh(roi_list):
    """Get dose-volume histogram curves from current plan.

    Parameters
    ----------
    roi_list : list of str
        Regions of interest to include in results.

    Returns
    -------
    dict
        Dose and volumes for given regions of interest.

    """
    dose = connect.get_current('Plan').TreatmentCourse.TotalDose
    max_dose = max([dose.GetDoseStatistic(RoiName=roi, DoseType='Max')
                    for roi in roi_list])
    dvh_dict = {'Dose': np.linspace(0, max_dose, 100)}
    for roi in roi_list:
        vals = dose.GetRelativeVolumeAtDoseValues(RoiName=roi,
                                                  DoseValues=dvh_dict['Dose'])
        dvh_dict[roi] = vals
    return dvh_dict
