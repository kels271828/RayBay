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

TODO: Add ways to compare results between plans (e.g., dvh plots, bars)

"""
import re

import numpy as np
import pandas as pd

import analyze


class RaybayResult:
    """RayStation treatment plan results.

    Attributes
    ----------
    patient : str
        Patient name.
    case : str
        Case name.
    plan : str
        Plan name.
    func_df : pandas.DataFrame
        Constituent function specifications.
    goal_df : pandas.DataFrame
        Clinical goal specificaions.
    roi_list : list of str
        Regions of interest included in clinical goals.
    norm : (str, float, float)
        Region of interest, dose, and volume used for normalization.
    utility : {'linear', 'linear_quadratic'}
        Shape of treatment plan utility function.
    solver : {'gp_minimize', 'forest_minimize', 'dummy_minimize'}
        Name of scikit-optimize solver used.
    time : float
        Total time in seconds for treatment plan optimization.
    flag_list : list
        RayStation exit statuses.
    opt_result : scipy.optimize.OptimizeResult
        Optimization results.
    goal_dict : dict
        Clinical goal results.
    dvh_dict : dict
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
    def __init__(self, patient, case, plan, funcs, norm, goals=None,
                 utility=None, solver=None):
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
        goals : str, optional
            Path to CSV with clinical goal specifications.
            If None, goals are assigned based on constituent functions.
        utility : {'linear', 'linear_quadratic'}, optional
            Shape of treatment plan utility function.
        solver : {'gp_minimize', 'forest_minimize', 'dummy_minimize'}, optional
            Name of scikit-optimize solver used.

        Returns
        -------
        RaybayResult
            RayStation treatment plan results.

        """
        self.patient = patient
        self.case = case
        self.plan = plan
        self.func_df = get_funcs(funcs)
        if isinstance(goals, str):
            self.goal_df = pd.read_csv(goals)
        else:
            self.goal_df = get_goals(self.funcs)
        self.roi_list = set(self.goal_df['Roi'])
        self.norm = norm
        self.utility = utility
        self.solver = solver
        self.time = 0.0
        self.flag_list = []
        self.opt_result = None
        self.goal_dict = {ii: [] for ii in range(len(self.goal_df))}
        self.dvh_dict = {}

    def get_utility(self, util_type=None):
        """Get treatment plan utility.

        Parameters
        ----------
        goal_df : pandas.DataFrame
            Clinical goal specifications.
        goal_dict : dict
            Clinical goal values.
        util_type : {'linear', 'linear_quadratic'}, optional
            Shape of treatment plan utility.

        Returns
        -------
        np.ndarray
            Vector of treatment plan utility values.

        """
        util_type = self.utility if util_type is None else util_type
        return utility(self.goal_df, self.goal_dict, util_type)

    def boxplot(self, data_type='goals', title=None, ax=None):
        """Visualize parameter and goal value ranges with a boxplot.

        Parameters
        ----------
        data_type : {'goals', 'pars'}, optional
            Type of boxplot to create.
        title : str, optional
            Figure title.
        ax : matplotlib.axes.Axes, optional
            Add the boxplot to the given axes.

        Returns
        -------
        None.

        """
        if data_type == 'pars':
            par_list = self.opt_result.x_iters
            analyze.boxplot(self.func_df, par_list, 'pars', title, ax)
        else:
            analyze.boxplot(self.goal_df, self.goal_dict, 'goals', title, ax)

    def corrplot(self, data_type='goals', title=None, size=500):
        """Visualize goal and parameter correlations with a heatmap.

        Modified from https://github.com/dylan-profiler/heatmaps.

        If data_type is 'pars', plots goals on the vertical axis and
        parameters on the horizontal axis. Otherwise plots goals on both
        vertical and horizontal axes.

        Parameters
        ----------
        data_type : {'goals', 'pars'}, optional
            Type of corrplot to create.
        title : str, optional
            Figure title.
        size : int, optional
            Size scale for boxes.

        Returns
        -------
        None.

        """
        if data_type == 'pars':
            analyze.corrplot(self.goal_df, self.goal_dict, self.func_df,
                             self.opt_result.x_iters, title=title, size=size)
        else:
            analyze.corrplot(self.goal_df, self.goal_dict, title=title,
                             size=size)

    def scatterplot(self, data_type='goals'):
        """Visualize goal and parameter relationships wiht scatterplots.

        If data_type is 'pars', plots goals on the vertical axis and
        parameters on the horizontal axis. Otherwise plots goals on both
        vertical and horizontal axes.

        Parameters
        ----------
        data_type : {'goals', 'pars'}, optional
            Type of scatterplot to create.

        Returns
        -------
        None.

        """
        if data_type == 'pars':
            analyze.scatterplot(self.goal_df, self.goal_dict, self.func_df,
                                self.opt_result.x_iters)
        else:
            analyze.scatterplot(self.goal_df, self.goal_dict)

    def dvhplot(self, roi_list=None):
        """Plot dose-volume histogram of solution.

        Parameters
        ----------
        roi_list : list of str, optional
            Regions of interest to include in figure.
            If None, all regions are included.

        Returns
        -------
        None.

        """
        roi_list = self.roi_list if roi_list is None else roi_list
        analyze.dvhplot(self.dvh_dict, roi_list)


def get_funcs(funcs):
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
    funcs : str
        Path to CSV with constituent function specifications.

    Returns
    -------
    pandas.DataFrame
        Constituent function specifications.

    """
    func_df = pd.read_csv(funcs).astype(object)
    for index, row in func_df.iterrows():
        for col in ['DoseLevel', 'PercentVolume', 'Weight']:
            # Tunable parameters are read in as strings '[min, max]',
            # so we need to convert them back to a list of floats.
            if isinstance(row[col], str):
                pars = [float(par) for par
                        in re.findall(r'\d+\.\d+|\d+', row[col])]
                func_df.loc[index, col] = pars if len(pars) > 1 else pars[0]
    return func_df


def get_goals(func_df):
    """Create clinical goals based on constituent functions.

    Clinical goals are specified by columns Roi, Type, GoalCriteria,
    AcceptanceLevel, and ParameterValue. Valid types include MinDose,
    AverageDose, MaxDose, MinDose, MinDvh, and MaxDvh.

    Parameters
    ----------
    func_df : pandas.DataFrame
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
    } for _, row in func_df.iterrows()]
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


def utility(goal_df, goal_dict, util_type):
    """Get treatment plan utility values.

    Parameters
    ----------
    goal_df : pandas.DataFrame
        Clinical goal specifications.
    goal_dict : dict
        Clinical goal values.
    util_type : {'linear', 'linear_quadratic'}
        Shape of treatment plan utility.

    Returns
    -------
    np.ndarray
        Vector of treatment plan utility values.
    """
    util_vec = np.zeros(len(goal_dict[0]))
    for ii in range(len(util_vec)):
        for index, row in goal_df.iterrows():
            util_vec[ii] += get_term(goal_dict[index][ii],
                                     row['AcceptanceLevel'], row['Type'],
                                     util_type)
    return util_vec


def get_term(value, level, goal_type, util_type):
    """Get treatment plan utility term value.

    Parameters
    ----------
    value : float
        Clinical goal value.
    level : float
        Clinical goal AcceptanceLevel.
    goal_type : str
        Clinical goal type (e.g., 'MaxDose')
    util_type : {'linear', 'linear_quadratic'}
        Shape of treatment plan utility term.

    Returns
    -------
    float
        Treatment plan utility term value.

    """
    if util_type not in ('linear', 'linear_quadratic'):
        raise ValueError(f'Invalid util_type: {util_type}')
    diff = 100*(value - level)/level
    if util_type == 'linear':
        return -diff if 'Max' in goal_type else diff
    else:
        if 'Max' in goal_type:
            return -diff if value <= level else -(diff + 1)*diff
        else:
            return diff if value >= level else -(diff - 1)*diff
