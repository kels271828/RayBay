import matplotlib.pyplot as plt
import numpy as np

import connect


def save_results(roi_names, fpath):
    """Save dose statistics and dvh curves from current plan.

    Parameters
    ----------
    roi_names : list
        Regions of interest to include in results.
    fpath : str
        File path to save results.

    """
    np.save(fpath + 'stats.npy', get_stats(roi_names))
    np.save(fpath + 'dvh.npy', get_dvh(roi_names))


def get_stats(roi_names):
    """Get dose statistics from current plan.

    Parameters
    ----------
    roi_names : list
        Regions of interest to include in results.

    Returns
    -------
    dict
        Dose statistics for given regions of interest.

    """
    stats = {}
    dose = connect.get_current('Plan').TreatmentCourse.TotalDose
    volumes = [0.99, 0.98, 0.95, 0.5, 0.05, 0.02, 0.01]
    volume_names = ['D99', 'D98', 'D95', 'D50', 'D5', 'D2', 'D1']
    for roi in roi_names:
        stats[roi] = {'Min': dose.GetDoseStatistic(RoiName=roi, DoseType='Min'),
                      'Average': dose.GetDoseStatistic(RoiName=roi,
                                                       DoseType='Average'),
                      'Max': dose.GetDoseStatistic(RoiName=roi, DoseType='Max')}
        doses = dose.GetDoseAtRelativeVolumes(RoiName=roi,
                                              RelativeVolumes=volumes)
        for ii in range(len(volumes)):
            stats[roi][volume_names[ii]] = doses[ii]
    return stats


def get_dvh(roi_names):
    """Get dvh curves from current plan.

    Parameters
    ----------
    roi_names : list
        Regions of interest to include in results.

    Returns
    -------
    dict
        Dose and volumes for given regions of interest.

    """
    dose = connect.get_current('Plan').TreatmentCourse.TotalDose
    max_dose = max([dose.GetDoseStatistic(RoiName=roi, DoseType='Max')
                    for roi in roi_names])
    dvh = {'Dose': np.linspace(0, max_dose, 100)}
    for roi in roi_names:
        dvh[roi] = dose.GetRelativeVolumeAtDoseValues(RoiName=roi,
                                                      DoseValues=dvh['Dose'])
    return dvh


def plot_dvh(roi_names):
    """Plot dvh curves from current plan.

    Parameters
    ----------
    roi_names : list
        Regions of interest to include in results.
        
    """
    dvh = get_dvh(roi_names)
    for roi in roi_names:
        plt.plot(dvh['Dose'], dvh[roi])
    plt.xlabel('Dose (cGy)')
    plt.ylabel('Volume (%)')
    plt.legend(roi_names)
