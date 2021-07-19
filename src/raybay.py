"""RayStation treatment planning with Bayesian optimization.

Constituent function parameters should be specified by columns Roi,
FunctionType, DoseLevel, PercentVolume, EudParameterA, and Weight. The
row index within the table should correspond to the constituent
function in the RayStation objective. Fixed parameters should be a
single value, tunable parameters should be a list containing the
minimum and maximum values, and irrelevant parameters can be left blank.

Clinical goals are specified by columns Roi, Type, GoalCriteria,
AcceptanceLevel, ParameterValue, Weight, and Shape. Valid types include
AverageDose, MinDose, MaxDose, MinDvh, and MaxDvh. Valid shapes include
linear and linear_quadratic.

Treatment plan utility function terms are specified by columms Weight
and Shape, and are included with the clinical goals. Valid shapes
include linear and linear_quadratic.

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
                 solver=None):
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
            self.goal_df = get_goals(self.func_df)
        self.roi_list = set(self.goal_df['Roi'])
        self.norm = norm
        self.solver = solver
        self.time = 0.0
        self.flag_list = []
        self.opt_result = None
        self.goal_dict = {ii: [] for ii in range(len(self.goal_df))}
        self.dvh_dict = {}

    def boxplot(self, data_type='goals', title=None):
        """Visualize parameter and goal value ranges with a boxplot.

        Parameters
        ----------
        data_type : {'goals', 'pars'}, optional
            Type of boxplot to create.
        flags : bool, optional
            If True, filter data by RayStation exit status.
        title : str, optional
            Figure title.

        Returns
        -------
        None.

        """
        if data_type == 'pars':
            par_list = self.opt_result.x_iters
            analyze.boxplot(self.func_df, par_list, 'pars', title)
        else:
            analyze.boxplot(self.goal_df, self.goal_dict, 'goals', title)

    def corrplot(self, data_type='goals', size=50, title=None):
        """Visualize goal and parameter correlations with a heatmap.

        Modified from https://github.com/dylan-profiler/heatmaps.

        If data_type is 'pars', plots goals on the vertical axis and
        parameters on the horizontal axis. Otherwise plots goals on both
        vertical and horizontal axes.

        Parameters
        ----------
        data_type : {'goals', 'pars'}, optional
            Type of corrplot to create.
        size : int, optional
            Size scale for boxes.
        title : str, optional
            Figure title.

        Returns
        -------
        None.

        """
        if data_type == 'pars':
            analyze.corrplot(self.goal_df, self.goal_dict, self.func_df,
                             self.opt_result.x_iters, size, title)
        else:
            analyze.corrplot(self.goal_df, self.goal_dict, size=size,
                             title=title)

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


class OptimizeResult:
    """Grid search parameter values.

    Because grid search doesn't have a utility function, we do not need
    to store opt_results as a scipy.optimize.OptimizeResult (with
    attributes like x, fun, func_vals, etc.). However, we want to be
    able to access x_iters similarly for plotting routines.

    Attributes
    ----------
    x_iters : list of lists
        Parameter evaluations for each iteration.

    """
    def __init__(self, x_iters):
        """Initialize instance of OptimizeResult

        Parameters
        ----------
        x_iters : list of lists
            Parameter evaluations for each iteration.

        Returns
        -------
        OptimizeResult
            Grid search parameter values.

        """
        self.x_iters = x_iters


def create_funcs(patient_path, case):
    """Format output of get_volumes.py into case-specific funcs.csv.

    If case == 'approved', CSV is initialized, but left blank. Run
    optimize.get_funcs() to get clinical constituent functions.
    Assumes all term weights are one. See get_dose_range() for
    additional tuneable parameter assumptions.

    Parameters
    ----------
    patient_path : str
        Path to patient folder.
    case : str
        Case name.

    Returns
    -------
    None.

    """
    if case == 'approved':
        func_df = pd.DataFrame(data={
            'Roi': [],
            'FunctionType': [],
            'DoseLevel': [],
            'PercentVolume': [],
            'EudParameterA': [],
            'Weight': []
        })
    else:
        roi_df = pd.read_csv(patient_path + 'goals.csv')
        n_roi = len(roi_df)
        func_df = pd.DataFrame(data={
            'Roi': roi_df['Roi'],
            'FunctionType': roi_df['Type'],
            'DoseLevel': roi_df['DoseLevel (cGy)'],
            'PercentVolume': roi_df['Volume (%)'],
            'EudParameterA': n_roi*[np.nan],
            'Weight': n_roi*[1]
        })
        func_df['PercentVolume'] = func_df.apply(add_zeros, axis=1)
        if case == 'bayes':
            func_df['DoseLevel'] = func_df.apply(get_dose_range, axis=1)
    func_df.sort_values(by='Roi', inplace=True)
    func_df.to_csv(patient_path + case + '/funcs.csv', index=False)


def add_zeros(row):
    """Replace missing PercentVolume values with zeros.

    Parameters
    ----------
    row : pandas.core.series.Series
        Row of func_df DataFrame.

    Returns
    -------
    float
        PercentVolume values.

    """
    if np.isnan(row['PercentVolume']):
        return 0
    return row['PercentVolume']


def get_dose_range(row):
    """Get tuneable parameter range.

    Assumes PTV D95 (MinDvh) is not tuneable.
    Assumes all other dose parameters (MaxDvh and MaxDose) are
    tuneable between a range of 0.25-1.0 of DoseLevel, except for
    PTV MaxDose, which uses (5*MaxDose - 4800)/4.

    Parameters
    ----------
    row : pandas.core.series.Series
        Row of func_df DataFrame.

    Returns
    -------
    float or str
        If parameter not tuneable, return DoseLevel.
        Otherwise, return lower and upper limits of DoseLevel.

    """
    if 'PTV' in row['Roi']:
        if row['FunctionType'] == 'MinDvh':
            return row['DoseLevel']
        min_dose = (row['DoseLevel'] + 3*4800)/4
        return f"[{min_dose}, {row['DoseLevel']}]"
    return f"[{row['DoseLevel']/4}, {row['DoseLevel']}]"


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


def create_goals(patient_path, case):
    """Format output of get_volumes.py into case-specific goals.csv.

    If case == 'bayes', goals.csv includes `Weight` and `Shape`
    columns with default values of 1 and 'linear'.

    Parameters
    ----------
    patient_path : str
        Path to patient folder.
    case : str
        Case name.

    Returns
    -------
    None.

    """
    roi_df = pd.read_csv(patient_path + 'goals.csv')
    goal_df = pd.DataFrame(data={
        'Roi': roi_df['Roi'],
        'Type': roi_df['Type'],
        'GoalCriteria': roi_df['GoalCriteria'],
        'AcceptanceLevel': roi_df['DoseLevel (cGy)'],
        'ParameterValue': roi_df['Volume (%)']
    })
    if case == 'bayes':
        goal_df['Weight'] = len(roi_df)*[1]
        goal_df['Shape'] = goal_df.apply(get_util_shape, axis=1)
    goal_df.sort_values(by='Roi', inplace=True)
    goal_df.to_csv(patient_path + case + '/goals.csv', index=False)


def get_util_shape(row):
    """Get utility term shape based on ROI.

    Parameters
    ----------
    row : pandas.core.series.Series
        Row of func_df DataFrame.

    Returns
    -------
    str
        If 'chest' or 'rib in row['Roi'], then return 'linear'.
        Otherwise, return 'linear_quadratic'.

    """
    if any([roi in row['Roi'].lower() for roi in ['chest', 'rib']]):
        return 'linear'
    return 'linear_quadratic'


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


def get_utility(goal_df, goal_dict, weights=None, shapes=None):
    """Get treatment plan utility values.

    Parameters
    ----------
    goal_df : pandas.DataFrame
        Clinical goal specifications.
    goal_dict : dict
        Clinical goal values.
    weights : list of float, optional
        Utility term weights. If None, uses Weight column in goal_df.
    shapes : list of str, optional
        Shape of utility terms ('linear' or 'linear_quadratic').
        If None, uses Shape column in goal_df.

    Returns
    -------
    np.ndarray
        Vector of treatment plan utility values.

    """
    if weights is None:
        weights = goal_df['Weight']
    if shapes is None:
        shapes = goal_df['Shape']
    n_util = len(goal_dict[0])
    util_vec = np.zeros(n_util)
    for ii in range(n_util):
        for index, row in goal_df.iterrows():
            util_vec[ii] += weights[index]*get_term(
                goal_dict[index][ii],
                row['AcceptanceLevel'],
                row['Type'], shapes[index])
    return util_vec


def get_term(value, level, goal_type, shape):
    """Get treatment plan utility term value.

    Parameters
    ----------
    value : float
        Clinical goal value.
    level : float
        Clinical goal AcceptanceLevel.
    goal_type : str
        Clinical goal type (e.g., 'MaxDose')
    shape : {'linear', 'linear_quadratic'}
        Shape of treatment plan utility term.

    Returns
    -------
    float
        Treatment plan utility term value.

    """
    if shape not in ('linear', 'linear_quadratic'):
        raise ValueError(f'Invalid shape: {shape}')
    diff = 100*(value - level)/level
    if shape == 'linear':
        return -diff if 'Max' in goal_type else diff
    if 'Max' in goal_type:
        return -diff if value <= level else -(diff + 1)*diff
    return diff if value >= level else -(diff - 1)*diff
