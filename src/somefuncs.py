"""Useful functions.

Note: Documentation has nicer ways to access some things, should update

What does ResetOptimization() do? Is it necessary to run each time?"""
import matplotlib.pyplot as plt
import numpy as np


def get_objective_terms(patient):
    """Get list of objective function term objects.
    
    Parameters
    ----------
    patient : connect.connect_cpython.PyScriptObject
        Current patient object.
    
    Returns
    -------
    list
        List of objective function term objects.

    """
    return patient.Cases[0].TreatmentPlans[1].PlanOptimizations[0].Objective.\
           ConstituentFunctions
           

def solve_oar_pars(patient, oar_term, par_type, par_vals, print_results=True):    
    """Solve treatment plan for various OAR parameters.
    
    Parameters
    ----------
    patient : connect.connect_cpython.PyScriptObject
        Current patient object.
    oar_term : TYPE
        OAR objective function term object.
    par_type : {'Weight', 'Dose', 'Volume'}
        OAR parameter type.
    par_vals : numpy.array
       Array of OAR parameter values.    
    print_results : bool, optional
        If True, print optimization results.
        
    Returns
    -------
    numpy.array
        Array of OAR average dose values.
    numpy.array
        Array of PTV D95 dose values.

    """
    oar_avg = np.zeros_like(par_vals)
    ptv_d95 = np.zeros_like(par_vals)
    for k in range(len(par_vals)):
        results = solve_oar_par(patient, oar_term, par_type, par_vals[k],
                                print_results=print_results)
        oar_avg[k] = results[0]
        ptv_d95[k] = results[1]
    return oar_avg, ptv_d95


def solve_oar_par(patient, oar_term, par_type, par_val,
                      print_results=True):
    """Solve treatment plan for OAR parameter.    

    Parameters
    ----------
    patient : connect.connect_cpython.PyScriptObject
        Current patient object.
    oar_term : TYPE
        OAR objective function term.
    par_type : {'Weight', 'Dose', 'Volume'}
        OAR parameter type.
    par_val : float
        OAR parameter value.    
    print_results : bool, optional
        If True, print optimization results.

    Returns
    -------
    float
        OAR average dose value.
    float
        PTV D95 dose value.

    """
    if print_results:
        print(f'{par_type}: {par_val:e}', end=', ')
    set_parameter(oar_term, par_type, par_val)
    run_optimization(patient)
    oar_avg = get_dose_stat(patient, get_roi_name(oar_term), 'Average')
    ptv_d95 = get_relative_dose(patient, 'PTV', 0.95)
    if print_results:
        print(f'{get_roi_name(oar_term)} Avg: {oar_avg:e}, PTV D95: {ptv_d95:e}')
    return oar_avg, ptv_d95


def set_parameter(oar_term, par_type, par_val):
    """Set OAR objective term parameter.
    
    Parameters
    ----------
    oar_term : connect.connect_cpython.PyScriptObject
        OAR objective function term.
    par_type : {'Weight', 'Dose', 'Volume'}
        OAR paramter type.
    par_val : float
        Oar parameter value.
        
    Returns
    -------
    None.
    
    """
    if par_type == 'Weight':
        oar_term.DoseFunctionParameters.Weight = par_val
    elif par_type == 'Dose':
        oar_term.DoseFunctionParameters.DoseLevel = par_val
    elif par_type == 'Volume':
        oar_term.DoseFunctionParameters.PercentVolume = par_val
    else:
        raise ValueError('Incorrect parameter type.')


def run_optimization(patient, reset=True):
    """Run optimization.

    Parameters
    ----------
    patient : connect.connect_cpython.PyScriptObject
        Current patient object.
    reset : bool, optional
        If True, reset optimization before running.

    Returns
    -------
    None.

    """
    if reset:
        reset_optimization(patient)
    patient.Cases[0].TreatmentPlans[1].PlanOptimizations[0].RunOptimization()


def reset_optimization(patient):
    """Reset optimization. 
    
    Parameters
    ----------
    patient : connect.connect_cpython.PyScriptObject
        Current patient object.

    Returns
    -------
    None.

    """
    patient.Cases[0].TreatmentPlans[1].PlanOptimizations[0].\
        ResetOptimization()


def get_roi_name(roi_term):
    """Get ROI name.
    
    Parameters
    ----------
    roi_term : connect.connect_cpython.PyScriptObject
        ROI objective function term object.

    Returns
    -------
    str
        ROI name.

    """
    return roi_term.ForRegionOfInterest.Name


def get_dose_stat(patient, roi_name, dose_type):
    """Get ROI dose statistic.    

    Parameters
    ----------
    patient : connect.connect_cpython.PyScriptObject
        Current patient object.
    roi_name : str
        ROI name.
    dose_type : {'Min', 'Average', 'Max'}
        Dose statistic type.

    Returns
    -------
    float
        ROI dose statistic.

    """
    return patient.Cases[0].TreatmentPlans[1].PlanOptimizations[0].\
        PlanningPhaseDose.GetDoseStatistic(RoiName=roi_name,
                                            DoseType=dose_type)
        
        
def get_relative_dose(patient, roi_name, volume):
    """Get ROI dose at relative volumes.

    Parameters
    ----------
    patient : connect.connect_cpython.PyScriptObject
        Current patient object.
    roi_name : str
        Region of interest name.
    vol_list : list
        List of relative volumes.

    Returns
    -------
    list
        List of ROI dose values.

    """
    return patient.Cases[0].TreatmentPlans[1].PlanOptimizations[0].\
        PlanningPhaseDose.GetDoseAtRelativeVolumes(RoiName=roi_name,
                                                  RelativeVolumes=[volume])[0]


def plot_results(par_vals, avg_vals, d95_vals, par_type, oar_name):
    """Plot optimization results.    

    Parameters
    ----------
    par_vals : numpy.array
        Array of OAR parameter values.
    avg_vals : numpy.array
        Array of OAR average dose values.
    d95_vals : numpy.array
        Array of PTV D95 dose values.
    par_type : {'Weight', 'Dose', 'Volume'}
        OAR parameter type.
    oar_name : str
        OAR name.

    Returns
    -------
    matplotlib.figure.Figure
        matplotlib figure object.
    matplotlib.axes.Axes
        matplotlib axes object.

    """
    fig, ax = plt.subplots(1, 3, figsize=(20, 5))
    ax[0].plot(par_vals, avg_vals, 'o')
    ax[0].set_xlabel(f'{oar_name} {par_type}')
    ax[0].set_ylabel(f'{oar_name} Average Dose')
    ax[1].plot(par_vals, d95_vals, 'o')
    ax[1].set_xlabel(f'{oar_name} {par_type}')
    ax[1].set_ylabel('PTV D95 Dose')
    ax[2].plot(avg_vals, d95_vals, 'o')
    ax[2].set_xlabel(f'{oar_name} Average Dose')
    ax[2].set_ylabel('PTV D95 Dose')
    return fig, ax
