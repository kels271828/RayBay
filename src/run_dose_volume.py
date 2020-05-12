"""Solve treatment plan for various OAR doses and volumes."""
import numpy as np

import connect
import somefuncs


# Get RayStation objects
patient = connect.get_current('Patient')
obj_terms = somefuncs.get_objective_terms(patient)
lung_term = obj_terms[11]

# Solve for different OAR dose and volume parameters
doses = np.arange(0.0, 3000.0, 500.0)
volumes = np.arange(0.0, 11.0, 2.0)
oar_avg = np.zeros((len(doses), len(volumes)))
ptv_d95 = np.zeros((len(doses), len(volumes)))
for ii in range(len(doses)):
    print(f'Dose: {doses[ii]:e}')
    somefuncs.set_parameter(lung_term, 'Dose', doses[ii])
    for jj in range(len(volumes)):
        print(f'    Volume: {volumes[jj]:e}', end=', ')
        somefuncs.set_parameter(lung_term, 'Volume', volumes[jj])
        try:
            somefuncs.run_optimization(patient)
            oar_avg[ii, jj] = somefuncs.get_dose_stat(patient, 'Lungs', 'Average')
            ptv_d95[ii, jj] = somefuncs.get_relative_dose(patient, 'PTV', 0.95)
            print(f'OAR Avg: {oar_avg[ii, jj]:e}, PTV D95: {ptv_d95[ii, jj]:e}')
        except Exception as e:
            oar_avg[ii, jj] = -1.0
            ptv_d95[ii, jj] = -1.0
            print('Something wrong: ' + str(e))

# Save results
np.save('doses.npy', doses)
np.save('volumes.npy', volumes)
np.save('oar_avg_dv.npy', oar_avg)
np.save('ptv_d95_dv.npy', ptv_d95)
