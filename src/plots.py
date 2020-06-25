"""Visualize sampled treatment plan results."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set(color_codes=True, font_scale=1.2)


def filter_flag(data, results, flag, keep=True):
    """Filter data by flag.

    Parameters
    ----------
    data : pandas.DataFrame
        Either sampled constituent function parameters
        or clinical goal results.
    results : pandas.DataFrame
        Clinical goal results.
    flag : {0, 1, 2}
        0 = success, 1 = normalization failed, 2 = optimization failed
    keep : bool, optional
        If True, keep values with given flag.
        If False, remove values with given flag.

    Returns
    -------
    pandas.DataFrame
        Data filtered by flag.

    """
    samples = results[results['Flag'] == flag]['Sample'].values
    return data[data['Sample'].isin(samples) == keep]


def boxplot(specs, values, data_type, title=None, ax=None):
    """Visualize parameter and goal value ranges with a boxplot.

    Parameters
    ----------
    specs : pandas.DataFrame
        Either constituent function specifications or
        clinical goal specifications.
    values : pandas.DataFrame
        Either sampled constituent function parameters or
        clinical goal results.
    data_type : {'pars', 'goals'}
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


def format_data(specs, values, data_type):
    """Format data and labels for boxplot and scatterplot.

    Parameters
    ----------
    specs : pandas.DataFrame
        Either constituent function specifications or
        clinical goal specifications.
    values : pandas.DataFrame
        Either sampled constituent function parameters or
        clinical goal results.
    data_type : {'pars', 'goals'}
        Type of data to format.

    Returns
    -------
    list
        Result values.
    list
        Result labels.

    """
    if data_type not in ('pars', 'goals'):
        raise ValueError(f'Invalid data_type: {data_type}')
    data, labels = [], []
    if data_type == 'pars':
        pars = get_tune_pars(specs)
        for _, row in pars.iterrows():
            data.append(values[values['Term'] == row['Term']][row['Par']])
            labels.append(row['Roi'] + ' ' + row['Par'])
    else:
        for index, row in specs.iterrows():
            data.append(values[index])
            labels.append(row['Roi'] + ' ' + row['Type'])
    return data, labels


def get_tune_pars(funcs):
    """Get tunable function parameters.

    The tunable function parameters are returned as a DataFrame with
    columns corresponding to each parameter: Term, Roi, and Par. The
    term column corresponds to the rows in the constituent function
    DataFrame.

    Parameters
    ----------
    funcs : pandas.DataFrame
        Constituent function specifications.

    Returns
    -------
    pandas.DataFrame
        Sampled tunable function parameters.

    """
    pars = []
    for idx, row in funcs.iterrows():
        for par in ['DoseLevel', 'PercentVolume', 'Weight']:
            if isinstance(row[par], list) or \
               (isinstance(row[par], str) and '[' in row[par]):
                pars.append({'Term': idx, 'Roi': row['Roi'], 'Par': par})
    return pd.DataFrame(data=pars, columns=['Term', 'Roi', 'Par'])


def corrplot(goals, results, funcs=None, pars=None, title=None):
    """Visualize goal and parameter correlations with a heatmap.

    Modified from https://github.com/dylan-profiler/heatmaps.

    If funcs and pars given, plots goals on the vertical axis and
    parameters on the horizontal axis. Otherwise plots goals on both
    vertical and horizontal axes.

    Parameters
    ----------
    goals : pandas.DataFrame
        Clinical goal specifications.
    results : pandas.DataFrame
        Clinical goal results.
    funcs : pandas.DataFrame, optional
        Constituent function specifications.
    pars : pandas.DataFrame, optional
        Sampled constituent function parameters.
    title : str, optional
        Figure title.

    Returns
    -------
    None.

    """
    # Format data
    ydata, ylabels = format_data(goals, results, 'goals')
    if funcs is None:
        xdata, xlabels = ydata, ylabels
    else:
        xdata, xlabels = format_data(funcs, pars, 'pars')
    xdata, xlabels = xdata[::-1], xlabels[::-1]

    # Plot boxes
    palette = sns.diverging_palette(20, 220, n=256)
    plot_grid = plt.GridSpec(1, 15)
    ax = plt.subplot(plot_grid[:, :-1])
    for ii in range(len(xdata)):
        for jj in range(len(ydata)):
            corr = np.corrcoef(xdata[ii], ydata[jj])[0, 1]
            ax.scatter(ii, jj, marker='s', s=1e3*abs(corr),
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


def scatterplot(goals, results, funcs=None, pars=None):
    """Visualize goal and parameter relationships with scatterplots.

    If funcs and pars given, plots goals on the vertical axis and
    parameters on the horizontal axis. Otherwise plots goals on both
    vertical and horizontal axes.

    Parameters
    ----------
    goals : pandas.DataFrame
        Clinical goal specifications.
    results : pandas.DataFrame
        Clinical goal results.
    funcs : pandas.DataFrame, optional
        Constituent function specifications.
    pars : pandas.DataFrame, optional
        Sampled constituent function parameters.

    Returns
    -------
    None.

    """
    ydata, ylabels = format_data(goals, results, 'goals')
    if funcs is None:
        xdata, xlabels = ydata, ylabels
    else:
        xdata, xlabels = format_data(funcs, pars, 'pars')
    for ii in range(len(ydata)):
        level = goals.iloc[ii]['AcceptanceLevel']
        fig, ax = plt.subplots(1, len(xdata), figsize=(25, 5))
        for jj in range(len(xdata)):
            ax[jj].plot(xdata[jj], ydata[ii], '.')
            ax[jj].plot([min(xdata[jj]), max(xdata[jj])], [level, level])
            ax[jj].set_xlabel(xlabels[jj])
            corr = f'Corr: {np.corrcoef(xdata[jj], ydata[ii])[0, 1]:.2f}'
            ax[jj].set_title(corr)
        ax[0].set_ylabel(ylabels[ii])
