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
doses = np.arange(1900.0, 2001.0, 100.0)
volumes = np.arange(0.0, 11.0, 1.0)
oar_avg = -np.ones((len(doses), len(volumes)))
ptv_max = -np.ones((len(doses), len(volumes)))
ptv_d95 = -np.ones((len(doses), len(volumes)))
dvh_doses = -np.ones((len(doses), len(volumes), 100))
dvh_oar = -np.ones((len(doses), len(volumes), 100))
dvh_ptv = -np.ones((len(doses), len(volumes), 100))
oar_avg_norm = -np.ones((len(doses), len(volumes)))
ptv_max_norm = -np.ones((len(doses), len(volumes)))
ptv_d95_norm = -np.ones((len(doses), len(volumes)))
dvh_doses_norm = -np.ones((len(doses), len(volumes), 100))
dvh_oar_norm = -np.ones((len(doses), len(volumes), 100))
dvh_ptv_norm = -np.ones((len(doses), len(volumes), 100))
fpath = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\results\\'
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
            
            # Store doses
            oar_avg[ii, jj] = plan_dose.GetDoseStatistic(RoiName='Lungs', DoseType='Average')
            ptv_max[ii, jj] = plan_dose.GetDoseStatistic(RoiName='PTV', DoseType='Max')
            ptv_d95[ii, jj] = plan_dose.GetDoseAtRelativeVolumes(RoiName='PTV', RelativeVolumes=[0.95])[0]

            # Store DVH
            max_oar = plan_dose.GetDoseStatistic(RoiName='Lungs', DoseType='Max')
            max_dose = max(ptv_max[ii, jj], max_oar)
            dvh_doses[ii, jj, :] = np.linspace(0, max_dose, num=100)
            dvh_ptv[ii, jj, :] = plan_dose.GetRelativeVolumeAtDoseValues(RoiName='PTV', DoseValues=dvh_doses[ii, jj, :])
            dvh_oar[ii, jj, :] = plan_dose.GetRelativeVolumeAtDoseValues(RoiName='Lungs', DoseValues=dvh_doses[ii, jj, :])
            
            print(f'OAR Avg: {oar_avg[ii, jj]:e}, PTV D95: {ptv_d95[ii, jj]:e}')
            
            # Store normalized doses
            beam_set.NormalizeToPrescription(RoiName='PTV', DoseValue=4800.0,
                                             DoseVolume=95.0, PrescriptionType='DoseAtVolume')
            oar_avg_norm[ii, jj] = plan_dose.GetDoseStatistic(RoiName='Lungs', DoseType='Average')
            ptv_max_norm[ii, jj] = plan_dose.GetDoseStatistic(RoiName='PTV', DoseType='Max')
            ptv_d95_norm[ii, jj] = plan_dose.GetDoseAtRelativeVolumes(RoiName='PTV', RelativeVolumes=[0.95])[0]
            
            # Store normalized DVH
            max_oar = plan_dose.GetDoseStatistic(RoiName='Lungs', DoseType='Max')
            max_dose = max(ptv_max_norm[ii, jj], max_oar)
            dvh_doses_norm[ii, jj, :] = np.linspace(0, max_dose, num=100)
            dvh_ptv_norm[ii, jj, :] = plan_dose.GetRelativeVolumeAtDoseValues(RoiName='PTV', DoseValues=dvh_doses_norm[ii, jj, :])
            dvh_oar_norm[ii, jj, :] = plan_dose.GetRelativeVolumeAtDoseValues(RoiName='Lungs', DoseValues=dvh_doses_norm[ii, jj, :])
        except Exception as e:
            print('Something wrong: ' + str(e))

    # Save results
    np.save(fpath + 'doses_5_28.npy', doses)
    np.save(fpath + 'volumes_5_28.npy', volumes)
    np.save(fpath + 'oar_avg_5_28.npy', oar_avg)
    np.save(fpath + 'ptv_max_5_28.npy', ptv_max)
    np.save(fpath + 'ptv_d95_5_28.npy', ptv_d95)
    np.save(fpath + 'dvh_doses_5_28.npy', dvh_doses)
    np.save(fpath + 'dvh_oar_5_28.npy', dvh_oar)
    np.save(fpath + 'dvh_ptv_5_28.npy', dvh_ptv)
    np.save(fpath + 'oar_avg_norm_5_28.npy', oar_avg_norm)
    np.save(fpath + 'ptv_max_norm_5_28.npy', ptv_max_norm)
    np.save(fpath + 'ptv_d95_norm_5_28.npy', ptv_d95_norm)
    np.save(fpath + 'dvh_doses_norm_5_28.npy', dvh_doses_norm)
    np.save(fpath + 'dvh_oar_norm_5_28.npy', dvh_oar_norm)
    np.save(fpath + 'dvh_ptv_norm_5_28.npy', dvh_ptv_norm)
