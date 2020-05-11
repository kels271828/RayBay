"""Solve treatment plan for various OAR volumes."""
import numpy as np

import connect
import somefuncs


# Get RayStation objects
patient = connect.get_current('Patient')
obj_terms = somefuncs.get_objective_terms(patient)
lung_term = obj_terms[11]

# Solve for different lung weights
lung_volumes = np.arange(0, 10.0, 2.0)
lung_vals, ptv_vals = somefuncs.solve_oar_pars(patient, lung_term, 'Volume',
                                               lung_volumes)

# Save results
np.save('lung_volumes.npy', lung_volumes)
np.save('lung_avgs_v.npy', lung_vals)
np.save('ptv_vals_v.npy', ptv_vals)

# Plot results
fig, ax = somefuncs.plot_results(lung_volumes, lung_vals, ptv_vals, 'Volume',
                                 'Lung')
fig.savefig('lung_volumes.png')
