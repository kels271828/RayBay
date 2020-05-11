"""Solve treatment plan for various OAR doses."""
import numpy as np

import connect
import somefuncs


# Get RayStation objects
patient = connect.get_current('Patient')
obj_terms = somefuncs.get_objective_terms(patient)
lung_term = obj_terms[11]

# Solve for different lung weights
lung_doses = np.arange(0, 3000.0, 500.0)
lung_vals, ptv_vals = somefuncs.solve_oar_pars(patient, lung_term, 'Dose',
                                               lung_doses)

# Save results
np.save('lung_doses.npy', lung_doses)
np.save('lung_avgs_d.npy', lung_vals)
np.save('ptv_vals_d.npy', ptv_vals)

# Plot results
fig, ax = somefuncs.plot_results(lung_doses, lung_vals, ptv_vals, 'Dose',
                                 'Lung')
fig.savefig('lung_doses.png')
