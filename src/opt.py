"""Solve treatment plan with Bayesian optimization."""
import numpy as np
import skopt

import connect
import sample


def gp_minimize(funcs, roi, dose, volume, goals=None, fpath='', random_state=0,
                n_calls=25, verbose=True):
    """Bayesian optimization using Gaussian Processes.

    Parameters
    ----------
    funcs : pandas.DataFrame or str
        Constituent function specifications or path to CSV file.
    roi : str
        Region of interest for normalization.
    dose : float
        Dose value for normalization.
    volume : float
        Dose volume for normalization.
    goals : pandas.DataFrame or str, optional
        Clinical goal specifications or path to CSV file.
        If None, goals are based on constituent functions.
    fpath : str, optional
        Path to save results.
        If not specified, results are saved to the current directory.
    random_state : int, optional
        Set random state to something other than None for reproducible results.
    n_calls : int, optional
        Number of calls to objective.
    verbose : bool, optional
        Control the verbosity.

    Returns
    -------
    scipy.optimize.OptimizeResult
        Optimization results.

    """
    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Define functions and goals
    funcs, goals, _, _, _, _ = sample.init_prob(funcs, goals)

    # Optimize
    obj = lambda pars: objective(plan, beam_set, funcs, goals, roi, dose,
                                 volume, pars)
    results = skopt.gp_minimize(obj, dimensions=get_dimensions(funcs),
                                random_state=random_state, n_calls=n_calls,
                                verbose=verbose)

    # Save results
    np.save(fpath + 'x.npy', results.x)
    np.save(fpath + 'x_iters.npy', results.x_iters)
    np.save(fpath + 'fun.npy', results.fun)
    np.save(fpath + 'func_vals.npy', results.func_vals)


def forest_minimize(funcs, roi, dose, volume, goals=None, fpath='',
                    random_state=0, n_calls=25, verbose=True):
    """Bayesian optimization using decision trees.

    Parameters
    ----------
    funcs : pandas.DataFrame or str
        Constituent function specifications or path to CSV file.
    roi : str
        Region of interest for normalization.
    dose : float
        Dose value for normalization.
    volume : float
        Dose volume for normalization.
    goals : pandas.DataFrame or str, optional
        Clinical goal specifications or path to CSV file.
        If None, goals are based on constituent functions.
    fpath : str, optional
        Path to save results.
        If not specified, results are saved to the current directory.
    random_state : int, optional
        Set random state to something other than None for reproducible results.
    n_calls : int, optional
        Number of calls to objective.
    verbose : bool, optional
        Control the verbosity.

    Returns
    -------
    scipy.optimize.OptimizeResult
        Optimization results.

    """
    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Define functions and goals
    funcs, goals, _, _, _, _ = sample.init_prob(funcs, goals)

    # Optimize
    obj = lambda pars: objective(plan, beam_set, funcs, goals, roi, dose,
                                 volume, pars)
    results = skopt.forest_minimize(obj, dimensions=get_dimensions(funcs),
                                    n_calls=n_calls, verbose=verbose,
                                    random_state=random_state)

    # Save results
    np.save(fpath + 'x.npy', results.x)
    np.save(fpath + 'x_iters.npy', results.x_iters)
    np.save(fpath + 'fun.npy', results.fun)
    np.save(fpath + 'func_vals.npy', results.func_vals)


def dummy_minimize(funcs, roi, dose, volume, goals=None, fpath='',
                   random_state=0, n_calls=25, verbose=True):
    """Random search by uniform sampling within given bounds.

    Parameters
    ----------
    funcs : pandas.DataFrame or str
        Constituent function specifications or path to CSV file.
    roi : str
        Region of interest for normalization.
    dose : float
        Dose value for normalization.
    volume : float
        Dose volume for normalization.
    goals : pandas.DataFrame or str, optional
        Clinical goal specifications or path to CSV file.
        If None, goals are based on constituent functions.
    fpath : str, optional
        Path to save results.
        If not specified, results are saved to the current directory.
    random_state : int, optional
        Set random state to something other than None for reproducible results.
    n_calls : int, optional
        Number of calls to objective.
    verbose : bool, optional
        Control the verbosity.

    Returns
    -------
    scipy.optimize.OptimizeResult
        Optimization results.

    """
    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Define functions and goals
    funcs, goals, _, _, _, _ = sample.init_prob(funcs, goals)

    # Optimize
    obj = lambda pars: objective(plan, beam_set, funcs, goals, roi, dose,
                                 volume, pars)
    results = skopt.dummy_minimize(obj, dimensions=get_dimensions(funcs),
                                   n_calls=n_calls, verbose=verbose,
                                   random_state=random_state)

    # Save results
    np.save(fpath + 'x.npy', results.x)
    np.save(fpath + 'x_iters.npy', results.x_iters)
    np.save(fpath + 'fun.npy', results.fun)
    np.save(fpath + 'func_vals.npy', results.func_vals)


def get_dimensions(funcs):
    """Get constituent function parameter dimensions.

    Parameters
    ----------
    funcs : pandas.DataFrame
        Constituent function specifications.

    Returns
    -------
    list
        Constituent function parameter dimensions.

    """
    dimensions = []
    for _, row in funcs.iterrows():
        for par in ['DoseLevel', 'PercentVolume', 'Weight']:
            if isinstance(row[par], list):
                dimensions.append(row[par])
    return dimensions


def objective(plan, beam_set, funcs, goals, roi, dose, volume, pars):
    """Objective function for Bayesian optimization.

    Parameters
    ----------
    plan : connect.connect_cpytyon.PyScriptObject
        Current treatment plan.
    beam_set : connect.connect_cpytyon.PyScriptObject
        Current beam set.
    funcs : pandas.DataFrame
        Constituent function specifications.
    goals : pandas.DataFrame
        Clinical goal specifications.
    roi : str
        Region of interest for normalization.
    dose : float
        Dose value for normalization.
    volume : float
        Dose volume for normalization.
    pars : list
        OAR dose level and percent volume parameters.

    Returns
    -------
    float
        Treatment plan score.

    """
    set_pars(plan, funcs, pars)
    flag = sample.calc_plan(plan, beam_set, roi, dose, volume)
    if flag == 0:
        return get_penalty(plan, goals)
    else:
        return 1e6


def set_pars(plan, funcs, pars):
    """Set objective function parameters.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    funcs : pandas.DataFrame
        Constituent function specifications.
    pars : list
        Constituent function parameters.

    Returns
    -------
    None.

    """
    count = 0
    const_funcs = plan.PlanOptimizations[0].Objective.ConstituentFunctions
    for index, row in funcs.iterrows():
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


def get_penalty(plan, goals):
    """

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    goals : pandas.DataFrame
        Clinical goal specifications.

    Returns
    -------
    float
        Plan penalty value.

    """
    penalty = 0
    results = get_results(plan, goals)
    for index, row in goals.iterrows():
        level = row['AcceptanceLevel']
        value = results[index]
        penalty += (value - level)/level
    return penalty


def get_results(plan, goals):
    """Get clinical goal results.

    Clinical goals are specified by columns Roi, Type, GoalCriteria,
    AcceptanceLevel, and ParameterValue. Only MinDose, AverageDose,
    MaxDose, MinDvh, and MaxDvh are evaluated. All other clinical goals
    are set to NaN.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    goals : pandas.DataFrame
        Clinical goal specifications.

    Returns
    -------
    dict
        Clinical goal results.

    """
    dose = plan.TreatmentCourse.TotalDose
    results = {}
    for index, row in goals.iterrows():
        results[index] = sample.get_goal_result(dose, row)
    return results
