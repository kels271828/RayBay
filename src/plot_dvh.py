"""Import and plot dose-volume histograms"""
import matplotlib.pyplot as plt
import numpy as np

import connect

# Get RayStation objects
patient = connect.get_current('Patient')
case = connect.get_current('Case')
plan = connect.get_current('Plan')
beam_set = connect.get_current('BeamSet')

# Get optimization objects
plan_opt = plan.PlanOptimizations[0]
plan_dose = plan_opt.PlanningPhaseDose
obj_terms = plan_opt.Objective.ConstituentFunctions
oar_term = obj_terms[11]

# Optimize plan
oar_term.DoseFunctionParameters.DoseLevel = 300.0
oar_term.DoseFunctionParameters.PercentVolume = 6.0
plan_opt.ResetOptimization()
plan_opt.RunOptimization()

# Print results
oar_avg = plan_dose.GetDoseStatistic(RoiName='Lungs',
                                     DoseType='Average')
ptv_d95 = plan_dose.GetDoseAtRelativeVolumes(RoiName='PTV',
                                            RelativeVolumes=[0.95])[0]
print(f'OAR Avg: {oar_avg:e}, PTV D95: {ptv_d95:e}')

# Get dvh curves
max_ptv = plan_dose.GetDoseStatistic(RoiName='PTV', DoseType='Max')
max_oar = plan_dose.GetDoseStatistic(RoiName='Lungs', DoseType='Max')
max_dose = max(max_ptv, max_oar)
doses = np.linspace(0, max_dose, num=100)
dvh_ptv = plan_dose.GetRelativeVolumeAtDoseValues(RoiName='PTV',
                                                  DoseValues=doses)
dvh_oar = plan_dose.GetRelativeVolumeAtDoseValues(RoiName='Lungs',
                                                  DoseValues=doses)

# Plot dvh curves
plt.plot(doses, dvh_ptv)
plt.plot(doses, dvh_oar)
plt.xlabel('Dose (cGy)')
plt.ylabel('Relative Volume (%)')
plt.legend(['PTV', 'Lungs'])
plt.title('Raw Data')

# Normalize results
beam_set.NormalizeToPrescription(RoiName='PTV', DoseValue=4800.0,
                                 DoseVolume=95.0,
                                 PrescriptionType='DoseAtVolume')

# Print results
oar_avg = plan_dose.GetDoseStatistic(RoiName='Lungs',
                                     DoseType='Average')
ptv_d95 = plan_dose.GetDoseAtRelativeVolumes(RoiName='PTV',
                                            RelativeVolumes=[0.95])[0]
print(f'OAR Avg: {oar_avg:e}, PTV D95: {ptv_d95:e}')

# Get dvh curves
max_ptv = plan_dose.GetDoseStatistic(RoiName='PTV', DoseType='Max')
max_oar = plan_dose.GetDoseStatistic(RoiName='Lungs', DoseType='Max')
max_dose = max(max_ptv, max_oar)
doses_norm = np.linspace(0, max_dose, num=100)
dvh_ptv_norm = plan_dose.GetRelativeVolumeAtDoseValues(RoiName='PTV',
                                                       DoseValues=doses_norm)
dvh_oar_norm = plan_dose.GetRelativeVolumeAtDoseValues(RoiName='Lungs',
                                                       DoseValues=doses_norm)

# Plot dvh curves
plt.figure()
plt.plot(doses_norm, dvh_ptv_norm)
plt.plot(doses_norm, dvh_oar_norm)
plt.xlabel('Dose (cGy)')
plt.ylabel('Relative Volume (%)')
plt.legend(['PTV', 'Lungs'])
plt.title('Normalized Data')
