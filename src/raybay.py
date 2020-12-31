"""RayStation treatment planning with Bayesian optimization.

Constituent function parameters should be specified by columns Roi,
FunctionType, DoseLevel, PercentVolume, EudParameterA, and Weight. The
row index within the table should correspond to the constituent
function in the RayStation objective. Fixed parameters should be a
single value, tunable parameters should be a list containing the
minimum and maximum values, and irrelevant parameters can be left blank.

Clinical goals are specified by columns Roi, Type, GoalCriteria,
AcceptanceLevel, and ParameterValue. Valid types include MinDose,
AverageDose, MaxDose, MinDose, MinDvh, and MaxDvh.

"""
import re

import numpy as np
import pandas as pd


class RaybayResult:
    """Optimized RayStation treatment plan results.

    Attributes
    ----------
    patient : str
        Patient name.
    case : str
        Case name.
    plan : str
        Plan name.
    funcs : pandas.DataFrame
        Constituent function specifications.
    goals : pandas.DataFrame
        Clinical goal specificaions.
    roi_list : list of str
        Regions of interest included in clinical goals.
    norm : (str, float, float)
        Region of interest, dose, and volume used for normalization.
    solver : {'gp_minimize', 'forest_minimize', 'dummy_minimize'}
        Name of scikit-optimize solver used.
    opt_result : scipy.optimize.OptimizeResult
        Optimization results.
    goal_result : dict
        Clinical goal results.
    dvh_result : dict
        Dose-volume histograms of solution.

    Note: The optimization results returned as an OptimizeResult object.
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

    def __init__(self, patient, case, plan, funcs, norm, solver, goals=None):
        """Initialise instance of RaybayResult.

        Parameters
        ----------
        patient : str
            Patient name.
        case : str
            Case name.
        plan : str
            Plan name.
        funcs : str
            Path to CSV with constituent function specifications.
        norm : (str, float, float)
            Region of interest, dose, and volume used for normalization.
        solver : {'gp_minimize', 'forest_minimize', 'dummy_minimize'}
            Name of scikit-optimize solver used.
        goals : str, optional
            Path to CSV with clinical goal specifications.
            If None, goals are assigned based on constituent functions.

        Returns
        -------
        RaybayResult
            RaybayResult object with plan information.

        """
        self.patient = patient
        self.case = case
        self.plan = plan
        self.funcs = get_funcs(funcs)
        self.norm = norm
        self.solver = solver
        if isinstance(goals, str):
            self.goals = pd.read_csv(goals)
        else:
            self.goals = get_goals(self.funcs)
        self.goal_result = {ii: [] for ii in range(len(self.goals))}
        self.roi_list = set(self.goals['Roi'])


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
        Path to CSV with constituent function specifications.

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
    AcceptanceLevel, and ParameterValue. Valid types include MinDose,
    AverageDose, MaxDose, MinDose, MinDvh, and MaxDvh.

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
