"""RayStation treatment planning with Bayesian optimization.

TODO:
* Additional information about format of funcs, goals, etc.
* Continue importing and formatting all of my raystation-facing functions
* add grid_search function (1D for now)
* add goal results dictionary as output
* add get_dvh function

Big question about results: can I modify a dictionary in place????

"""
import re

import numpy as np
import pandas as pd
import skopt

import connect


def plan_opt(funcs, norm, goals=None, solver='gp_minimize', n_calls=25,
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
    funcs : pandas.DataFrame or str
        Constituent function specifications or path to CSV file.
    norm : (str, float, float)
        Region of interest, dose, and volume used for normalization.
    goals : pandas.DataFrame or str, optional
        Clinical goal specifications or path to CSV file.
        If None, goals are assigned based on constituent functions.
    solver : {'gp_minimize', 'forest_minimize', 'dummy_minimize'}
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
    scipy.optimize.OptimizeResult
        The optimization results returned as an OptimizeResult object.
        Important attributes are:

        - `x` [list]: location of the minimum.
        - `fun` [float]: function value at the minimum.
        - `models`: surrogate models used for each iteration.
        - `x_iters` [list of lists]: location of function evaluate for
           each iteration.
        - `func_vals` [array]: function value for each iteration.
        - `space` [Space]: the optimization space.
        - `specs` [dict]: the call specifications.
        - `rng` [RandomState instance]: State of the random state at the
          end of minimization.

        For more details about the OptimizeResult object, refer to
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html

    """
    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Define functions and goals
    if isinstance(funcs, str):
        funcs = get_funcs(funcs)
    if goals is None:
        goals = get_goals(funcs)
    elif isinstance(goals, str):
        goals = pd.read_csv(goals)

    # Optimize
    obj = lambda pars: objective(plan, beam_set, funcs, goals, norm, pars)
    if solver == 'forest_minimize':
        return skopt.forest_minimize(obj, dimensions=get_dimensions(funcs),
                                     n_calls=n_calls,
                                     n_initial_points=n_initial_points,
                                     random_state=random_state,
                                     verbose=verbose)
    elif solver == 'dummy_minimize':
        return skopt.dummy_minimize(obj, dimensions=get_dimensions(funcs),
                                    n_calls=n_calls,
                                    n_initial_points=n_initial_points,
                                    random_state=random_state,
                                    verbose=verbose)
    else:
        return skopt.gp_minimize(obj, dimensions=get_dimensions(funcs),
                                 n_calls=n_calls,
                                 n_initial_points=n_initial_points,
                                 random_state=random_state, verbose=verbose)


def get_funcs(fpath):
    """Load constituent functions from CSV file.

    Constituent function parameters should be specified by columns Roi,
    FunctionType, DoseLevel, PercentVolume, EudParameterA, and Weight.
    The row index within the table should correspond to the constituent
    function in the RayStation objective. Fixed parameters should be a
    a single value, tunable parameters should be a list containing the
    minimum and maximum values, and irrelevant parameters can be left
    blank.

    Parameters
    ----------
    fpath : str
        Path to CSV file.

    Returns
    -------
    pandas.DataFrame
        Constituent function specifications.

    """
    funcs = pd.read_csv(fpath).astype(object)
    for index, row in funcs.iterrows():
        for col in ['DoseLevel', 'PercentVolume', 'Weight']:
            # Tunable parameters are read in as strings '[min, max]',
            # so we need to convert them back to a list of floats.
            if isinstance(row[col], str):
                pars = [float(par) for par
                        in re.findall(r'\d+\.\d+|\d+', row[col])]
                funcs.loc[index, col] = pars if len(pars) > 1 else pars[0]
    return funcs


def get_goals(funcs):
    """Create clinical goals based on constituent functions.

    Clinical goals are specified by columns Roi, Type, GoalCriteria,
    AcceptanceLevel, and ParameterValue.

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
        'GoalCriteria': 'AtMost'
                        if 'Max' in row['FunctionType'] else 'AtLeast',
        'AcceptanceLevel': get_bound(row['DoseLevel'], row['FunctionType']),
        'ParameterValue': row['EudParameterA']
                          if 'Eud' in row['FunctionType'] else
                          get_bound(row['PercentVolume'], row['FunctionType'])
    } for _, row in funcs.iterrows()]
    columns = ['Roi', 'Type', 'GoalCriteria', 'AcceptanceLevel',
               'ParameterValue']
    return pd.DataFrame(data=data, columns=columns)


def get_bound(par, func_type):
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


def objective(plan, beam_set, funcs, goals, norm, pars):
    """Objective function for hyperparameter optimization.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    beam_set : connect.connect_cpython.PyScriptObject
        Current beam set.
    funcs : pandas.DataFrame
        Constituent function specifications.
    goals : pandas.DataFrame
        Clinical goal specifications.
    norm : (str, float, float)
        Region of interest, dose, and volume used for normalization.
    pars : list
        Constituent function parameters.

    Returns
    -------
    float
        Treatment plan score.

    """
    set_pars(plan, funcs, pars)
    flag = calc_plan(plan, beam_set, norm)
    print(f'Flag: {flag}')
    return get_score(plan, goals, flag)


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
        beam_set.NormalizeToPrescription(RoiName=norm[0], DoseValue=norm[1],
                                         DoseVolume=norm[2],
                                         PrescriptionType='DoseAtVolume')
        return 0
    except:
        return 1


def get_score(plan, goals, flag):
    """Calculate treatment plan score.

    The treatment plan score is a linear combination of the relative
    difference between goal values and goal results. Returns 1e6 if
    RayStation optimization or normalization failed.

    Parameters
    ----------
    plan : connect.connect_cpython.PyScriptObject
        Current treatment plan.
    goals : pandas.DataFrame
        Clinical goal specifications.
    flag : int
        RayStation exit status.

    Returns
    -------
    float
        Treatment plan score.

    """
    if flag > 0:
        return 1e6
    score = 0
    results = get_results(plan, goals)
    for index, row in goals.iterrows():
        level = row['AcceptanceLevel']
        value = results[index]
        sign = 1 if 'Most' in row['GoalCriteria'] else -1
        score += sign*(value - level)/level
    return score


def get_results(plan, goals):
    """Get clinical goal results.

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
