"""RayStation treatment planning with Bayesian optimization.

TODO:
* Add 1D grid_search function

Note: My virtual environments in Citrix have old version of skopt, so
kwarg n_initial_points reverted back to n_random_starts

"""
import re

import numpy as np
import skopt

import connect
import raybay


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
    funcs : str
        Path to CSV with constituent function specifications.
    norm : (str, float, float)
        Region of interest, dose, and volume used for normalization.
    goals : pandas.DataFrame or str, optional
        Path to CSV with clinical goal specifications.
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
    raybay.RaybayResult

    """
    # Get RayStation objects
    patient = connect.get_current('Patient')
    case = connect.get_current('Case')
    plan = connect.get_current('Plan')
    beam_set = connect.get_current('BeamSet')

    # Initialize result object
    result = raybay.RaybayResult(patient.Name, case.CaseName, plan.Name, funcs,
                                 norm, solver, goals)

    # Optimize
    def obj(pars):
        return objective(plan, beam_set, result.funcs, result.goals, norm,
                         result.goal_result, pars)
    if solver == 'forest_minimize':
        results = skopt.forest_minimize(obj, dimensions=get_dims(result.funcs),
                                        n_calls=n_calls,
                                        n_random_starts=n_initial_points,
                                        random_state=random_state,
                                        verbose=verbose)
    elif solver == 'dummy_minimize':
        results = skopt.dummy_minimize(obj, dimensions=get_dims(result.funcs),
                                       n_calls=n_calls,
                                       n_random_starts=n_initial_points,
                                       random_state=random_state,
                                       verbose=verbose)
    else:
        results = skopt.gp_minimize(obj, dimensions=get_dims(result.funcs),
                                    n_calls=n_calls,
                                    n_random_starts=n_initial_points,
                                    random_state=random_state, verbose=verbose)
    result.opt_result = results

    # Get optimal dose-volume histogram
    set_pars(plan, result.funcs, result.opt_result.x)
    calc_plan(plan, beam_set, norm)
    result.dvh_result = get_dvh(result.roi_list)

    return result


def objective(plan, beam_set, funcs, goals, norm, goal_result, pars):
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
    goal_result : dict
        Clinical goal results.
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
    return get_score(plan, goals, flag, goal_result)


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


def get_score(plan, goals, flag, goal_result):
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
    goal_result : dict
        Clinical goal results.

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
        goal_result[index].append(value)
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


def get_dims(funcs):
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
