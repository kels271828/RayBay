"""Solve for different doses"""
import matplotlib.pyplot as plt
import numpy as np

import connect

# Get patient
Patient = connect.get_current('Patient')

# Get objective function terms
ConstituentFunctions = Patient.Cases[0].TreatmentPlans[1].\
    PlanOptimizations[0].Objective.ConstituentFunctions

# Get specific terms
PTV_max = ConstituentFunctions[8]
PTV_min = ConstituentFunctions[9]
Lung_dvh = ConstituentFunctions[11]

# Solve for different lung weights
print('Solving for different doses')
dose_vec = np.arange(0, 3000.0, 500.0)
lung_avg = np.zeros_like(dose_vec)
ptv_d95 = np.zeros_like(dose_vec)
for k in range(len(dose_vec)):
    print('Lung dose:  %e' % dose_vec[k], end=', ')
    Lung_dvh.DoseFunctionParameters.DoseLevel = dose_vec[k]
    Patient.Cases[0].TreatmentPlans[1].PlanOptimizations[0].ResetOptimization()
    Patient.Cases[0].TreatmentPlans[1].PlanOptimizations[0].RunOptimization()
    lung_avg[k] = Patient.Cases[0].TreatmentPlans[1].PlanOptimizations[0].\
        PlanningPhaseDose.GetDoseStatistic(RoiName='Lungs',DoseType='Average')
    print('Lung avg: %e' % lung_avg[k], end=', ')
    ptv_d95[k] = Patient.Cases[0].TreatmentPlans[1].PlanOptimizations[0].\
        PlanningPhaseDose.\
        GetDoseAtRelativeVolumes(RoiName='PTV', RelativeVolumes=[0.95])
    print('PTV D95: %e' % ptv_d95[k])

# Plot results
fig, ax = plt.subplots(1, 3, figsize=(15, 5))
ax[0].plot(dose_vec, lung_avg, 'o')
ax[0].set_xlabel('Lung Weight')
ax[0].set_ylabel('Lung Average')
ax[1].plot(dose_vec, ptv_d95, 'o')
ax[1].set_xlabel('Lung Weight')
ax[1].set_ylabel('PTV D95')
ax[2].plot(lung_avg, ptv_d95, 'o')
ax[2].set_xlabel('Lung Average')
ax[2].set_ylabel('PTV D95')
