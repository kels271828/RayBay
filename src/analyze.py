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
import scipy.cluster.hierarchy as sch
import seaborn as sns
sns.set(color_codes=True, font_scale=1.2)


def boxplot(specs, values, data_type, title=None):
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

    Returns
    -------
    None.

    """
    data, labels = format_data(specs, values, data_type)
    fig, ax = plt.subplots(1, 1)
    ax.boxplot(data)
    ax.set_xticklabels(labels, rotation=90)
    if data_type == 'pars':
        ax.set_ylabel('Parameter Values')
    else:
        ax.set_ylabel('Goal Values')
    ax.set_title(title)


def corrplot(goal_df, goal_dict, func_df=None, par_list=None, size=50,
             cluster=True, title=None):
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
    size : int, optional
        Size scale for boxes.
    cluster : bool, optional
        If True, use hierarchical clustering to order rows and columns.
    title : str, optional
        Figure title.

    Returns
    -------
    None.

    """
    # Format data
    ydata, ylabels = format_data(goal_df, goal_dict, 'goals')
    if cluster:
        ydata, ylabels = cluster_data(goal_dict, ydata, ylabels, 'goals')
    if func_df is None:
        xdata, xlabels = ydata, ylabels
    else:
        xdata, xlabels = format_data(func_df, par_list, 'pars')
        if cluster:
            xdata, xlabels = cluster_data(goal_dict, xdata, xlabels, 'pars')
    ydata, ylabels = ydata[::-1], ylabels[::-1]  # reorder

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
            corr = np.corrcoef(xdata[jj], ydata[ii])[0, 1]
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


def goalplot(goal_df, res_dict, percent=True):
    """Visualize comparison of goal results as table.

    Parameters
    ----------
    goal_df : pandas.DataFrame
        Clinical goal specifications.
    res_dict : {str: raybay.RaybayResult}
        RayStation treatment plan names and results.
    percent : bool, optional
        If True, prints goal results as percentages, otherwise prints
        goal values.

    Returns
    -------
    None.

    """
    val_list = []
    roi_dict = {}
    for index, row in goal_df.iterrows():
        goal_vals = []
        level = row['AcceptanceLevel']
        roi_dict[f"{index} {row['Roi']} {row['Type']}"] = index
        if not percent:
            goal_vals.append(level)
        for res in res_dict:
            if res.lower() in ('approved', 'default'):
                idx_sol = 0
            else:
                idx_sol = np.argmin(res_dict[res].opt_result.func_vals)
            val = res_dict[res].goal_dict[index][idx_sol]
            if percent:
                goal_vals.append(100*(val - level)/level)
            else:
                goal_vals.append(val)
        val_list.append(goal_vals)
    if percent:
        columns = res_dict.keys()
    else:
        columns = ['AcceptanceLevel'] + list(res_dict.keys())
    val_df = pd.DataFrame(columns=columns, data=val_list, index=roi_dict)
    ncols = len(res_dict) if percent else 1.5*(len(res_dict) + 1)
    fig, ax = plt.subplots(1, figsize=(ncols, len(goal_df)))
    sns.heatmap(val_df, cmap=sns.diverging_palette(220, 20, n=256), center=0,
                annot=True, fmt=".2f", cbar_kws={'label': 'Goal Value'}, ax=ax)


def cluster_data(goal_dict, data, labels, data_type):
    """Get indices for hierarchical goal value clusters.

    Modified from https://github.com/TheLoneNut/CorrelationMatrixClustering/blob/master/CorrelationMatrixClustering.ipynb

    Parameters
    ----------
    goal_dict : dict
        Clinical goal results.
    data : list
        Result values.
    labels : list
        Result labels.
    data_type : {'goals', 'pars'}
        Type of data to cluster.

    Returns
    -------
    list
        Clustered result values.
    list
        Clustered result labels.

    """
    # Cluster goal results
    corr = pd.DataFrame(goal_dict).corr().values
    dist = sch.distance.pdist(corr)
    link = sch.linkage(dist, method='complete')
    idx_goals = np.argsort(sch.fcluster(link, 0.5*dist.max(), 'distance'))

    # Sort data
    if data_type == 'goals':
        data_sorted = np.array(data)[idx_goals]
        labels_sorted = [labels[idx] for idx in idx_goals]
    else:
        idx_pars = np.array([int(label.split()[0]) for label in labels])
        data_sorted = []
        labels_sorted = []
        for idx in idx_goals:
            if idx in idx_pars:
                ii = np.where(idx_pars == idx)[0][0]
                data_sorted.append(data[ii])
                labels_sorted.append(labels[ii])
    return data_sorted, labels_sorted


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
    flag : int, optional
        RayStation exit status.

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
            labels.append(f"{row['Term']} {row['Roi']} {row['Par']}")
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
