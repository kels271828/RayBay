"""Solve treatment plan for various OAR weights."""
import numpy as np

import connect
import somefuncs


# Get RayStation objects
patient = connect.get_current('Patient')
obj_terms = somefuncs.get_objective_terms(patient)
lung_term = obj_terms[11]

# Solve for different lung weights
lung_weights = 10.0**np.arange(0, 5)
lung_vals, ptv_vals = somefuncs.solve_oar_pars(patient, lung_term, 'Weight',
                                               lung_weights)

# Save results
np.save('lung_weights.npy', lung_weights)
np.save('lung_avgs_w.npy', lung_vals)
np.save('ptv_vals_w.npy', ptv_vals)

# Plot results
fig, ax = somefuncs.plot_results(lung_weights, lung_vals, ptv_vals, 'Weight',
                                 'Lung')
fig.savefig('lung_weights.png')
