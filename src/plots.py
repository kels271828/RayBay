"""Visualize sampled treatment plan results."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set(color_codes=True, font_scale=1.2)


def boxplot(specs, values, data_type, title=None):
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

    Returns
    -------
    None.

    """
    data, labels = _format_data(specs, values, data_type)
    plt.boxplot(data)
    plt.xticks(ticks=np.arange(1, len(labels) + 1), labels=labels, rotation=90)
    if data_type == 'pars':
        plt.ylabel('Parameter Values')
    else:
        plt.ylabel('Goal Vaues')
    if title is not None:
        plt.title(title)


def _format_data(specs, values, data_type):
    """Format data and labels for boxplot and scatterplot."""
    if data_type not in ('pars', 'goals'):
        raise ValueError(f'Invalid data_type: {data_type}')
    data, labels = [], []
    if data_type == 'pars':
        pars = _get_tune_pars(specs)
        for _, row in pars.iterrows():
            data.append(values[values['Term'] == row['Term']][row['Par']])
            labels.append(row['Roi'] + ' ' + row['Par'])
    else:
        for index, row in specs.iterrows():
            data.append(values[index])
            labels.append(row['Roi'] + ' ' + row['Type'])
    return data, labels


def _get_tune_pars(funcs):
    """Get tunable parameters."""
    pars = []
    for idx, row in funcs.iterrows():
        for par in ['DoseLevel', 'PercentVolume', 'Weight']:
            if isinstance(row[par], str) and '[' in row[par]:
                pars.append({'Term': idx, 'Roi': row['Roi'], 'Par': par})
    return pd.DataFrame(data=pars, columns=['Term', 'Roi', 'Par'])


def _filter(results, data, flag):
    """Filter data by flag."""
    samples = results[results['Flag'] == flag]['Sample'].values
    return data[data['Sample'].isin(samples)]


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
