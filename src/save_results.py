"""Save results for specific ROIs."""
import numpy as np

import connect

# ROIs
roi_names = ['SpinalCanal', 'Lungs', 'Lung_L', 'Lung_R', 'Heart',
             'Chestwall_L', 'Rib', 'PTV']

# Get RayStation objects
patient = connect.get_current('Patient')
case = connect.get_current('Case')
plan = connect.get_current('Plan')
beam_set = connect.get_current('BeamSet')

# Get optimization objects
plan_opt = plan.PlanOptimizations[0]
plan_dose = plan_opt.PlanningPhaseDose

# Get dose statistics
stats = {}
max_dose = 0
volumes = [0.99, 0.98, 0.95, 0.5, 0.05, 0.02, 0.01]
volume_names = ['D99', 'D98', 'D95', 'D50', 'D5', 'D2', 'D1']
for roi in roi_names:
    stats[roi] = {
        'Min': plan_dose.GetDoseStatistic(RoiName=roi, DoseType='Min'),
        'Average': plan_dose.GetDoseStatistic(RoiName=roi, DoseType='Average'),
        'Max': plan_dose.GetDoseStatistic(RoiName=roi, DoseType='Max')
    }
    max_dose = max(max_dose, stats[roi]['Max'])
    doses = plan_dose.GetDoseAtRelativeVolumes(RoiName=roi,
                                               RelativeVolumes=volumes)
    for ii in range(len(volumes)):
        stats[roi][volume_names[ii]] = doses[ii]

# Get dvh curves
dvh = {'Doses': np.linspace(0, max_dose, 100)}
for roi in roi_names:
    dvh[roi] = plan_dose.GetRelativeVolumeAtDoseValues(RoiName=roi,
                                                       DoseValues=dvh['Doses'])

# Save result
fpath = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\results\\'
np.save(fpath + 'stats.npy', stats)
np.save(fpath + 'dvh.npy', dvh)
