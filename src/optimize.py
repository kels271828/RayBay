"""RayStation treatment planning with Bayesian optimization."""
import pickle
import re
from time import time

import numpy as np
import pandas as pd
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


def grid_search(funcs, norm, goals=None, n_points=25, weight=False):
    """1D grid search for RayStation treatment planning.

    Parameters
    ----------
    funcs : str
        Path to CSV with constituent function specifications.
    norm : (str, float, float)
        Region of interest, dose, and volume used for normalization.
    goals : pandas.DataFrame or str, optional
        Path to CSV with clinical goal specifications.
        If None, goals are assigned based on constituent functions.
    n_points : int, optional
        Number of treatment plans to evaluate.
    weight : bool, optional
        If True, uses logspacing for grid values.

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
                                 norm, goals)
    dims = get_dims(result.func_df)[0]
    if weight:
        pars = np.geomspace(dims[0], dims[1], n_points)
    else:
        pars = np.linspace(dims[0], dims[1], n_points)

    # Evaluate treatment plans
    start_time = time()
    for par in pars:
        set_pars(plan, result.func_df, [par])
        flag = calc_plan(plan, beam_set, norm)
        result.flag_list.append(flag)
        print(f'Par: {par}, Flag: {flag}, Time: {time() - start_time:.4f}',
              flush=True)
        results = get_results(plan, result.goal_df)
        scale = get_scale(result.goal_df, result.norm, results) \
            if flag == 1 else 1.0
        for index, _ in result.goal_df.iterrows():
            value = scale*results[index]
            result.goal_dict[index].append(value)
        with open(funcs[:-9] + 'goal_dict.pkl', 'wb') as fp:
            pickle.dump(result.goal_dict, fp)
    result.time = time() - start_time

    # Save parameter values
    x_iters = [[par] for par in pars]
    result.opt_result = raybay.OptimizeResult(x_iters)

    return result


def grid_search2(funcs, norm, goals=None, n_points=[25, 25]):
    """2D grid search for RayStation treatment planning.

    Parameters
    ----------
    funcs : str
        Path to CSV with constituent function specifications.
    norm : (str, float, float)
        Region of interest, dose, and volume used for normalization.
    goals : pandas.DataFrame or str, optional
        Path to CSV with clinical goal specifications.
        If None, goals are assigned based on constituent functions.
    n_points : list of int, optional
        Number of treatment plans to evaluate for each dimension.

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
                                 norm, goals)
    dims = get_dims(result.func_df)[:2]
    pars = [np.linspace(dims[ii][0], dims[ii][1], n_points[ii])
            for ii in range(2)]
    x_iters = []

    # Evaluate treatment plans
    start_time = time()
    for ii in pars[0]:
        for jj in pars[1]:
            x_iters.append([ii, jj])
            set_pars(plan, result.func_df, [ii, jj])
            flag = calc_plan(plan, beam_set, norm)
            result.flag_list.append(flag)
            time_iter = time() - start_time
            print(f'Pars: {ii, jj}, Flag: {flag}, Time: {time_iter:.4f}',
                  flush=True)
            results = get_results(plan, result.goal_df)
            scale = get_scale(result.goal_df, result.norm, results) \
                if flag == 1 else 1.0
            for index, _ in result.goal_df.iterrows():
                value = scale*results[index]
                result.goal_dict[index].append(value)
            with open(funcs[:-9] + 'goal_dict.pkl', 'wb') as fp:
                pickle.dump(result.goal_dict, fp)
    result.time = time() - start_time

    # Save parameter values
    result.opt_result = raybay.OptimizeResult(x_iters)

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
    return get_utility(plan, result.goal_df, result.norm, flag,
                       result.goal_dict, repo_path)


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


def get_utility(plan, goal_df, norm, flag, goal_dict, repo_path):
    """Calculate treatment plan utility.

    Returns 1e6 if RayStation optimization failed.
    Negative of utility computed in raybay module (since skopt solves
    a minimization problem).

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
        Treatment plan utility.

    """
    if flag == 2:
        return 1e6
    results = get_results(plan, goal_df)
    scale = get_scale(goal_df, norm, results) if flag == 1 else 1.0
    utility = 0
    for index, row in goal_df.iterrows():
        value = scale*results[index]
        goal_dict[index].append(value)
        term = raybay.get_term(value, row['AcceptanceLevel'], row['Type'],
                               row['Shape'])
        utility += -row['Weight']*term
    with open(repo_path + 'goal_dict.pkl', 'wb') as fp:
        pickle.dump(goal_dict, fp)
    return utility


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
    if 'Dvh' in goal['Type']:
        volume = 0.01*goal['ParameterValue']
        return dose.GetDoseAtRelativeVolumes(RoiName=goal['Roi'],
                                             RelativeVolumes=[volume])[0]
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


def get_volumes(patient_path):
    """Get ROI names and volumes.

    Initalize patient goal spreadsheet with `Roi` and `Volume (cm^3)`
    columns. Saved as a CSV file in the provided directory.

    Parameters
    ----------
    patient_path : str
        Path to patient folder.

    Returns
    -------
    None.

    """
    # Get RayStation objects
    case = connect.get_current('Case')
    roi_geometries = case.PatientModel.StructureSets[0].RoiGeometries

    # Get names and volumes
    roi_names = []
    roi_volumes = []
    for roi in roi_geometries:
        roi_name = roi.OfRoi.Name
        try:
            roi_volume = roi.GetRoiVolume()
        except:
            roi_volume = 'Error'
        roi_names.append(roi_name)
        roi_volumes.append(roi_volume)

    # Save results
    n_roi = len(roi_geometries)
    roi_df = pd.DataFrame(data={
        'Roi': roi_names,
        'RoiVolume (cm^3)': roi_volumes,
        'Type': n_roi*[np.nan],
        'GoalCriteria': n_roi*[np.nan],
        'DoseLevel (cGy)': n_roi*[np.nan],
        'Volume (cm^3)': n_roi*[np.nan],
        'Volume (%)': n_roi*[np.nan]
    })
    roi_df.to_csv(patient_path + 'goals.csv', index=False)


def get_funcs(plan):
    """Get clinical constituent functions from plan.

    Currently only able to handle MinDose, MaxDose, MinDvh, MaxDvh,
    and DoseFall-Off function types. Does not extract EudParameter A
    values from plan.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.


    Returns
    -------
    pandas.DataFrame
        Constituent function specifications.

    """
    func_df = pd.DataFrame(data={
        'Roi': [],
        'FunctionType': [],
        'DoseLevel': [],
        'PercentVolume': [],
        'EudParameterA': [],
        'Weight': []
    })
    const_funcs = plan.PlanOptimizations[0].Objective.ConstituentFunctions
    for func in const_funcs:
        func_pars = func.DoseFunctionParameters
        try:
            func_df = func_df.append({
                'Roi': func.ForRegionOfInterest.Name,
                'FunctionType': func_pars.FunctionType,
                'DoseLevel': func_pars.DoseLevel,
                'PercentVolume': func_pars.PercentVolume,
                'Weight': func_pars.Weight
            }, ignore_index=True)
        except:
            high_dose = func_pars.HighDoseLevel
            low_dose = func_pars.LowDoseLevel
            dose_dist = func_pars.LowDoseDistance
            func_type = f"Dose Fall-Off [H]{high_dose} cGy "
            func_type += f"[L]{low_dose} cGy, "
            func_type += f"Low dose distance {dose_dist} cm"
            func_df = func_df.append({
                'Roi': func.ForRegionOfInterest.Name,
                'FunctionType': func_type,
                'Weight': func_pars.Weight
            }, ignore_index=True)
    return func_df


def add_funcs(plan, func_df):
    """Add constituent function terms to plan.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    func_df : pandas.DataFrame
        Constituent function specifications.

    Returns
    -------
    None.

    """
    plan_opt = plan.PlanOptimizations[0]
    plan_opt.ClearConstituentFunctions()
    for _, row in func_df.iterrows():
        plan_opt.AddOptimizationFunction(FunctionType=row['FunctionType'],
                                         RoiName=row['Roi'])
