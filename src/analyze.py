"""Visualize sampled treatment plan results.

TODO:
- add wrappers in raybay!
- what other functions would be useful?
- plot dvh (how many differnet plans/ROIS?)

"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set(color_codes=True, font_scale=1.2)


def boxplot(specs, values, data_type, title=None, ax=None):
    """Visualize parameter and goal value ranges with a boxplot.

    Parameters
    ----------
    specs : pandas.DataFrame
        Either constituent function specifications or
        clinical goal specifications.
    values : dict or list
        Either clinical goal results or sampled function parameters.
    data_type : {'goals', 'pars'}
        Type of boxplot to create.
    title: str, optional
        Figure title.
    ax : matploblib.axes.Axes, optional
        Add the boxplot to the given axes.

    Returns
    -------
    None.

    """
    data, labels = format_data(specs, values, data_type)
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    ax.boxplot(data)
    ax.set_xticklabels(labels, rotation=90)
    if data_type == 'pars':
        ax.set_ylabel('Parameter Values')
    else:
        ax.set_ylabel('Goal Vaues')
    ax.set_title(title)


def corrplot(goal_df, goal_dict, func_df=None, par_list=None, title=None,
             size=500):
    """Visualize goal and parameter correlations with a heatmap.

    Modified from https://github.com/dylan-profiler/heatmaps.

    If funcs and pars given, plots goals on the vertical axis and
    parameters on the horizontal axis. Otherwise plots goals on both
    vertical and horizontal axes.

    Parameters
    ----------
    goal_df : pandas.DataFrame
        Clinical goal specifications.
    goal_dict : dict
        Clinical goal results.
    func_df : pandas.DataFrame, optional
        Constituent function specifications.
    par_list : list, optional
        Sampled constituent function parameters.
    title : str, optional
        Figure title.
    size : int, optional
        Size scale for boxes.

    Returns
    -------
    None.

    """
    # Format data
    ydata, ylabels = format_data(goal_df, goal_dict, 'goals')
    if func_df is None:
        xdata, xlabels = ydata, ylabels
    else:
        xdata, xlabels = format_data(func_df, par_list, 'pars')
    xdata, xlabels = xdata[::-1], xlabels[::-1]

    # Plot boxes
    palette = sns.diverging_palette(20, 220, n=256)
    plot_grid = plt.GridSpec(1, 15)
    ax = plt.subplot(plot_grid[:, :-1])
    for ii in range(len(xdata)):
        for jj in range(len(ydata)):
            corr = np.corrcoef(xdata[ii], ydata[jj])[0, 1]
            ax.scatter(ii, jj, marker='s', s=size*abs(corr),
                       c=[palette[int(255/2*(corr + 1))]])

    # Initialize tick labels
    ax.set_xticks(range(len(xdata)))
    ax.set_xticklabels(xlabels, rotation=90)
    ax.set_yticks(range(len(ydata)))
    ax.set_yticklabels(ylabels)
    ax.set_title(title)

    # Adjust grid lines relative to boxes
    ax.grid(False, 'major')
    ax.grid(True, 'minor')
    ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
    ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)
    ax.set_xlim([-0.5, len(xlabels) - 0.5])
    ax.set_ylim([-0.5, len(ylabels) - 0.5])

    # Legend
    ax = plt.subplot(plot_grid[:, -1])
    bar_y = np.linspace(-1, 1, len(palette))
    bar_h = bar_y[1] - bar_y[0]
    ax.barh(y=bar_y, width=[1]*len(palette), height=bar_h, color=palette,
            linewidth=0)
    ax.set_facecolor('w')
    ax.set_xticks([])
    ax.set_yticks([-1, 0, 1])
    ax.yaxis.tick_right()


def scatterplot(goal_df, goal_dict, func_df=None, par_list=None):
    """Visualize goal and parameter relationships with scatterplots.

    If funcs and pars given, plots goals on the vertical axis and
    parameters on the horizontal axis. Otherwise plots goals on both
    vertical and horizontal axes.

    Parameters
    ----------
    goal_df : pandas.DataFrame
        Clinical goal specifications.
    goal_dict : dict
        Clinical goal results.
    func_df : pandas.DataFrame, optional
        Constituent function specifications.
    par_list : list, optional
        Sampled constituent function parameters.

    Returns
    -------
    None.

    """
    ydata, ylabels = format_data(goal_df, goal_dict, 'goals')
    if func_df is None:
        xdata, xlabels = ydata, ylabels
    else:
        xdata, xlabels = format_data(func_df, par_list, 'pars')
    for ii in range(len(ydata)):
        level = goal_df.iloc[ii]['AcceptanceLevel']
        fig, ax = plt.subplots(1, len(xdata), figsize=(25, 5))
        for jj in range(len(xdata)):
            ax[jj].plot(xdata[jj], ydata[ii], '.')
            ax[jj].plot([min(xdata[jj]), max(xdata[jj])], [level, level])
            ax[jj].set_xlabel(xlabels[jj])
            corr = f'Corr: {np.corrcoef(xdata[jj], ydata[ii])[0, 1]:.2f}'
            ax[jj].set_title(corr)
        ax[0].set_ylabel(ylabels[ii])


def format_data(specs, values, data_type):
    """Format data and labels for boxplot and scatterplot.

    Parameters
    ----------
    specs : pandas.DataFrame
        Either constituent function specifications or
        clinical goal specifications.
    values : dict or list
        Either clinical goal results or sampled function parameters.
    data_type : {'goals', 'pars'}
        Type of data to format.

    Returns
    -------
    list
        Result values.
    list
        Result labels.

    """
    if data_type not in ('goals', 'pars'):
        raise ValueError(f'Invalid data_type: {data_type}')
    data, labels = [], []
    if data_type == 'goals':
        for index, row in specs.iterrows():
            data.append(values[index])
            labels.append(row['Roi'] + ' ' + row['Type'])
    else:
        pars = get_pars(specs)
        for index, row in pars.iterrows():
            data.append([value[index] for value in values])
            labels.append(row['Roi'] + ' ' + row['Par'])
    return data, labels


def get_pars(func_df):
    """Get tunable function parameters.

    The tunable function parameters are returned as a DataFrame with
    columns corresponding to each parameter: Term, Roi, and Par. The
    term column corresponds to the rows in the constituent function
    DataFrame.

    Parameters
    ----------
    func_df : pandas.DataFrame
        Constituent function specifications.

    Returns
    -------
    pandas.DataFrame
        Sampled tunable function parameters.

    """
    pars = []
    for idx, row in func_df.iterrows():
        for par in ['DoseLevel', 'PercentVolume', 'Weight']:
            # Tunable parameters are read in as strings '[min, max]'
            # when funcs loaded from CSV rather than RaybayResult.
            if isinstance(row[par], list) or \
               (isinstance(row[par], str) and '[' in row[par]):
                pars.append({'Term': idx, 'Roi': row['Roi'], 'Par': par})
    return pd.DataFrame(data=pars, columns=['Term', 'Roi', 'Par'])
