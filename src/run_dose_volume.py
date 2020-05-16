"""Solve treatment plan for various OAR doses and volumes."""
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
oar_term = obj_terms[11] # Lungs

# Solve for different OAR dose and volume parameters
doses = np.arange(0.0, 2001.0, 100.0)
volumes = np.arange(0.0, 11.0, 1.0)
oar_avg = -np.ones((len(doses), len(volumes)))
ptv_d95 = -np.ones((len(doses), len(volumes)))
oar_avg_norm = -np.ones((len(doses), len(volumes)))
ptv_d95_norm = -np.ones((len(doses), len(volumes)))
for ii in range(len(doses)):
    print(f'Dose: {doses[ii]:e}')
    oar_term.DoseFunctionParameters.DoseLevel = doses[ii]
    for jj in range(len(volumes)):
        print(f'    Volume: {volumes[jj]:e}', end=', ')
        oar_term.DoseFunctionParameters.PercentVolume = volumes[jj]
        try:
            # Run optimization
            plan_opt.ResetOptimization()
            plan_opt.RunOptimization()
            
            # Save results
            oar_avg[ii, jj] = plan_dose.GetDoseStatistic(RoiName='Lungs',
                                                          DoseType='Average')
            ptv_d95[ii, jj] = plan_dose.GetDoseAtRelativeVolumes(RoiName='PTV',
                                                                 RelativeVolumes=[0.95])[0]
            
            # Save normalized results
            beam_set.NormalizeToPrescription(RoiName='PTV', DoseValue=4800.0,
                                             DoseVolume=95.0,
                                             PrescriptionType='DoseAtVolume')
            oar_avg_norm[ii, jj] = plan_dose.GetDoseStatistic(RoiName='Lungs',
                                                              DoseType='Average')
            ptv_d95_norm[ii, jj] = plan_dose.GetDoseAtRelativeVolumes(RoiName='PTV',
                                                                      RelativeVolumes=[0.95])[0]
            print(f'OAR Avg: {oar_avg[ii, jj]:e}, PTV D95: {ptv_d95[ii, jj]:e}')
        except Exception as e:
            print('Something wrong: ' + str(e))

# Save results
fpath = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\results\\'
np.save(fpath + 'doses_5_15.npy', doses)
np.save(fpath + 'volumes_5_15.npy', volumes)
np.save(fpath + 'oar_avg_5_15.npy', oar_avg)
np.save(fpath + 'ptv_d95_5_15.npy', ptv_d95)
