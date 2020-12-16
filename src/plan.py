"""RayStation treatment planning with Bayesian optimization.

TODO:
* Additional information about format of funcs, goals, etc.
* Continue importing and formatting all of my raystation-facing functions

"""
import re
import pickle

import numpy as np
import pandas as pd
import skopt

import connect


def plan_opt(funcs, norm, goals=None, solver='gp_minimize', fpath='',
             n_calls=25, random_state=None, n_initial_points=10,
             verbose=True):
    """Hyperparameter optimization for RayStation treatment planning.

    Hyperparameter optimization for RayStation treatment planning using
    the following solvers from scikit-optimize:

        - `gp_minimize`: Bayesian optimization using Gaussian Processes
        - `forest_minimize`: Sequential optimization using decision
           trees
        - `dummy_minimize`: Random search by uniform sampling within
           the given bounds

    For more details related to scikit-optimize, refer to
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
    solver : {'gp_minimize', 'forest_minimize', 'random_minimize'}
        Name of scikit-optimize solver to use.
    fpath : str, optional
        Path to output directory.
        If not specified, results are saved in the current directory.
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
        - `rng` [RandomState instance]: State of the random state at
          at the end of minimization.

        For more details related to the OptimizeResult object, refer to
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html

    """
    # Get RayStation objects
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Define functions and goals
    if isinstance(funcs, str):
        funcs = load_funcs(funcs)
    if goals is None:
        goals = init_goals(funcs)
    elif isinstance(goals, str):
        goals = pd.read_csv(goals)

    # Optimize
    # maybe create actual objective
    # option for different solvers
    obj = lambda pars: objective(plan, beam_set, funcs, goals, roi, dose,
                                 volume, pars)
    results = skopt.gp_minimize(obj, dimensions=get_dimensions(funcs),
                                random_state=random_state, n_calls=n_calls,
                                verbose=verbose)

    # Save results
    # why don't i just save the OptimizeResult?
    np.save(fpath + 'x.npy', results.x)
    np.save(fpath + 'x_iters.npy', results.x_iters)
    np.save(fpath + 'fun.npy', results.fun)
    np.save(fpath + 'func_vals.npy', results.func_vals)
    pickle.dump(results.models, open(fpath + 'models.npy', 'wb'))


def load_funcs(fpath):
    """Load constituent functions from CSV file.

    Constituent function parameters should be specified by columns
    Roi, FunctionType, DoseLevel, PercentVolume, EudParameterA, and
    Weight. The row index within the table should correspond to the
    constituent function in the RayStation objective. Fixed parameters
    should be a single value, tunable parameters should be a list
    containing the minimum and maximum values, and irrelevant
    parameters can be left blank.

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
    pandas.DataFrame
        Clinical goal specifications.

    """
    data = [{
        'Roi': row['Roi'],
        'Type': row['FunctionType'],
        'GoalCriteria': 'AtMost'
                        if 'Max' in row['FunctionType'] else 'AtLeast',
        'AcceptanceLevel': get_par_bound(row['DoseLevel'],
                                         row['FunctionType']),
        'ParameterValue': row['EudParameterA']
                          if 'Eud' in row['FunctionType'] else
                          get_par_bound(row['PercentVolume'],
                                        row['FunctionType'])
    } for _, row in funcs.iterrows()]
    columns = ['Roi', 'Type', 'GoalCriteria', 'AcceptanceLevel',
               'ParameterValue']
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
