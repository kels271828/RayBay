"""Solve for various OAR weights"""
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
print('Solving for different weights')
weight_vec = 10.0**np.arange(0, 5)
lung_avg = np.zeros_like(weight_vec)
ptv_d95 = np.zeros_like(weight_vec)
for k in range(len(weight_vec)):
    print('Lung weight:  %e' % weight_vec[k], end=', ')
    Lung_dvh.DoseFunctionParameters.Weight = weight_vec[k]
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
ax[0].plot(weight_vec, lung_avg, 'o')
ax[0].set_xlabel('Lung Weight')
ax[0].set_ylabel('Lung Average')
ax[1].plot(weight_vec, ptv_d95, 'o')
ax[1].set_xlabel('Lung Weight')
ax[1].set_ylabel('PTV D95')
ax[2].plot(lung_avg, ptv_d95, 'o')
ax[2].set_xlabel('Lung Average')
ax[2].set_ylabel('PTV D95')
