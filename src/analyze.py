"""Visualize sampled treatment plan results.

TODO:
- plot dvh
- compare dvh
- print result percent change
- compare result percent change

"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set(color_codes=True, font_scale=1.2)


def boxplot(specs, values, data_type, flag_list=None, title=None):
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
    flag_list : list, optional
        RayStation exit statuses.
    title: str, optional
        Figure title.

    Returns
    -------
    None.

    """
    data, labels = format_data(specs, values, data_type, flag_list)
    fig, ax = plt.subplots(1, len(data), squeeze=False,
                           figsize=(6.4*len(data), 4.8))
    for flag in range(len(data)):
        ax[0, flag].boxplot(data[flag])
        ax[0, flag].set_xticklabels(labels, rotation=90)
        if flag == 0:
            if data_type == 'pars':
                ax[0, flag].set_ylabel('Parameter Values')
            else:
                ax[0, flag].set_ylabel('Goal Values')
        if flag_list is None:
            ax[0, flag].set_title(title)
        else:
            ax[0, flag].set_title(f'Flag: {flag}')


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


def scatterplot(goal_df, goal_dict, func_df=None, par_list=None,
                flag_list=None):
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
    flag_list : list, optional
        RayStation exit statuses.

    Returns
    -------
    None.

    """
    ydata, ylabels = format_data(goal_df, goal_dict, 'goals', flag_list)
    if func_df is None:
        xdata, xlabels = ydata, ylabels
    else:
        xdata, xlabels = format_data(func_df, par_list, 'pars', flag_list)
    for ii in range(len(ydata[0])):
        level = goal_df.iloc[ii]['AcceptanceLevel']
        fig, ax = plt.subplots(1, len(xdata[0]), figsize=(25, 5))
        for jj in range(len(xdata[0])):
            xmin = []
            xmax = []
            for flag in range(len(xdata)):
                ax[jj].plot(xdata[flag][jj], ydata[flag][ii], '.')
                xmin.append(min(xdata[flag][jj]))
                xmax.append(max(xdata[flag][jj]))
            ax[jj].plot([min(xmin), max(xmax)], [level, level])
            ax[jj].set_xlabel(xlabels[jj])
            corr = np.corrcoef(xdata[flag][jj], ydata[flag][ii])[0, 1]
            ax[jj].set_title(f'Corr: {corr:.2f}')
        ax[0].set_ylabel(ylabels[ii])


def dvhplot(dvh_dict, roi_list):
    """Plot dose-volume histogram of solution.

    Parameters
    ----------
    dvh_dict : dict
        Dose-volume histograms of solution.
    roi_list : list of str
        Regions of interest to include in figure.

    Returns
    -------
    None.

    """
    for roi in roi_list:
        plt.plot(dvh_dict['Dose'], dvh_dict[roi])
    plt.xlabel('Dose (cGy)')
    plt.ylabel('Volume (%)')
    plt.legend(roi_list, bbox_to_anchor=(1, 1))


def format_data(specs, values, data_type, flag_list=None):
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
    flag_list : list, optional
        RayStation exit statuses.

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
            labels.append(f"{index} {row['Roi']} {row['Type']}")
    else:
        pars = get_pars(specs)
        for index, row in pars.iterrows():
            data.append([value[index] for value in values])
            labels.append(row['Roi'] + ' ' + row['Par'])
    return filter_flags(data, flag_list), labels


def filter_flags(data, flag_list):
    """Filter figure data by RayStation exit status.

    Parameters
    ----------
    data : list
        Result values.
    flag_list : list
        RayStation exit statuses.

    Returns : list
        Result values filtered by RayStation exit status.

    """
    if flag_list is None:
        return [data]
    data_flag = []
    for flag in set(flag_list):
        data_flag.append([[row[ii] for ii in range(len(flag_list))
                           if flag_list[ii] == flag] for row in data])
    return data_flag


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
