"""Visualize sampled treatment plan results."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_pars(funcs, pars, results, flag=0):
    """Plot parameter values as boxplots.

    Parameters
    ----------
    funcs : pandas.DataFrame
        Constituent function specifications.
    pars : pandas.DataFrame
        Sampled constituent function parameters.
    results : pandas.DataFrame
        Clinical goal results.
    flag : int, optional
        0 = success, 1 = normalization failed, 2 = optimization failed

    Returns
    -------
    None.

    """
    data, labels = _get_par_data(funcs, _filter(results, pars, flag))
    plt.boxplot(data)
    plt.xticks(ticks=np.arange(1, len(labels) + 1), labels=labels, rotation=90)
    plt.ylabel('Parameter Value')


def _filter(results, data, flag):
    """Filter data by flag."""
    samples = results[results['Flag'] == flag]['Sample'].values
    return data[data['Sample'].isin(samples)]


def _get_par_data(funcs, pars):
    """Get data and labels for parameter boxplot."""
    tune_pars = _get_tune_pars(funcs)
    data = [pars[pars['Term'] == row['Term']][row['Parameter']]
            for _, row in tune_pars.iterrows()]
    labels = [row['Roi'] + ' ' + row['Parameter']
              for _, row in tune_pars.iterrows()]
    return data, labels


def _get_tune_pars(funcs):
    """Get tunable parameters."""
    tune_pars = pd.DataFrame(columns=['Term', 'Parameter'])
    for index, row in funcs.iterrows():
        for col in ['DoseLevel', 'PercentVolume', 'EudParameterA', 'Weight']:
            if isinstance(row[col], str) and '[' in row[col]:
                tune_pars = tune_pars.append({'Term': index, 'Roi': row['Roi'],
                                              'Parameter': col},
                                             ignore_index=True)
    return tune_pars


def plot_goals(goals, results, flag=0):
    """Plot goal values as boxplots.

    Parameters
    ----------
    goals : pandas.DataFrame
        Clinical goal specifications.
    results : pandas.DataFrame
        Clinical goal results.
    flag : int, optional
        0 = success, 1 = normalization failed

    Returns
    -------
    None.

    """
    data, labels = _get_goal_data(goals, results[results['Flag'] == flag])
    plt.boxplot(data)
    plt.xticks(ticks=np.arange(1, len(labels) + 1), labels=labels, rotation=90)
    plt.ylabel('Goal Value')


def _get_goal_data(goals, results):
    """Get data and labels for goal boxplot."""
    data = [results[ii] for ii in range(len(goals))]
    labels = [row['Roi'] + ' ' + row['Type']
              for _, row in goals.iterrows()]
    return data, labels


def plot_goals_pars(funcs, pars, goals, results, flag=0):
    """Plot goals v. parameters as scatter plots.

    Parameters
    ----------
    funcs : pandas.DataFrame
        Constituent function specifications.
    pars : pandas.DataFrame
        Sampled constituent function parameters.
    goals : pandas.DataFrame
        Clinical goal specifications.
    results : pandas.DataFrame
        Clinical goal results
    flag : int, optional
        0 = success, 1 = normalization failed
    Returns
    -------
    None.

    """
    par_data, par_labels = _get_par_data(funcs, _filter(results, pars, flag))
    goal_data, goal_labels = _get_goal_data(goals,
                                            results[results['Flag'] == flag])
    for ii in range(len(goal_data)):
        fig, ax = plt.subplots(1, len(par_data), figsize=(25, 5))
        level = goals.iloc[ii]['AcceptanceLevel']
        for jj in range(len(par_data)):
            ax[jj].plot(par_data[jj], goal_data[ii], '.')
            ax[jj].plot([min(par_data[jj]), max(par_data[jj])], [level, level])
            ax[jj].set_xlabel(par_labels[jj])
            corr = f'Corr: {np.corrcoef(par_data[jj], goal_data[ii])[0, 1]:.2}'
            ax[jj].set_title(corr)
        ax[0].set_ylabel(goal_labels[ii])


def plot_corr(data, title='', x_labels=None, y_labels=None):
    """Plot correlations by color and size.

    Modified from https://github.com/dylan-profiler/heatmaps.

    Parameters
    ----------
    data : pandas.DataFrame
        Results to include in plot.
    title : str, optional
        Plot title.
    x_labels : list, optional
        Columns to appear on horizontal axis.
        If none, all columns are included.
    y_labels : list, optional
        Columns to appear on vertical axis.
        If none, all columns are included.

    Returns
    -------
    None.

    """
    # Format and filter data
    corr = data.corr()
    if x_labels is not None:
        corr = corr[x_labels]
    if y_labels is not None:
        corr = corr.loc[y_labels]
    corr = pd.melt(corr.reset_index(), id_vars='index').replace(np.nan, 0)
    corr.columns = ['x', 'y', 'z']
    x = corr['x']
    y = corr['y']
    x_names = [t for t in sorted(set([v for v in x]))]
    x_to_num = {p[1]: p[0] for p in enumerate(x_names)}
    y_names = [t for t in sorted(set([v for v in y]))]
    y_to_num = {p[1]: p[0] for p in enumerate(y_names)}

    palette = sns.diverging_palette(20, 220, n=256)

    def _value_to_size(val):
        val_position = 0.99*val + 0.01
        val_position = min(max(val_position, 0), 1)
        return 100*val_position

    def _value_to_color(val):
        val_position = float((val + 1))/2
        val_position = min(max(val_position, 0), 1)
        ind = int(val_position*(len(palette) - 1))
        return palette[ind]

    # Plot squares
    plot_grid = plt.GridSpec(1, 15, hspace=0.2, wspace=0.1)
    ax = plt.subplot(plot_grid[:, :-1])
    ax.scatter(
        x=[x_to_num[v] for v in x],
        y=[y_to_num[v] for v in y],
        marker='s',
        s=[_value_to_size(v) for v in corr['z'].abs()],
        c=[_value_to_color(v) for v in corr['z']]
    )

    # Annotations
    ax.set_xticks([v for k, v in x_to_num.items()])
    ax.set_xticklabels([k for k in x_to_num], rotation=45,
                       horizontalalignment='right')
    ax.set_yticks([v for k, v in y_to_num.items()])
    ax.set_yticklabels([k for k in y_to_num])

    ax.grid(False, 'major')
    ax.grid(True, 'minor')
    ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
    ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)

    ax.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5])
    ax.set_ylim([-0.5, max([v for v in y_to_num.values()]) + 0.5])
    ax.set_facecolor('#F1F1F1')

    ax.set_title(title)

    # Legend
    ax = plt.subplot(plot_grid[:, -1])

    col_x = [0]*len(palette)
    bar_y = np.linspace(-1, 1, len(palette))

    bar_height = bar_y[1] - bar_y[0]
    ax.barh(
        y=bar_y,
        width=[5]*len(palette),
        left=col_x,
        height=bar_height,
        color=palette,
        linewidth=0
    )
    ax.set_xlim(1, 2)
    ax.grid(False)
    ax.set_facecolor('white')
    ax.set_xticks([])
    ax.set_yticks(np.linspace(min(bar_y), max(bar_y), 3))
    ax.yaxis.tick_right()
