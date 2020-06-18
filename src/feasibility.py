"""
Figure out feasible parameters!
"""
import matplotlib.pyplot as plt
import numpy as np

import connect
import bayes2

# Get RayStation objects
plan = connect.get_current('Plan')
beam_set = connect.get_current('BeamSet')

# Set ROIs
roi_names = ['PTV 4/7/20', 'Lungs', 'SpinalCord (Thorax)',
             'Esophagus', 'Heart']
max_vals = [6271, 2000, 5000, 6930, 3500]

# Look for percent difference with lowest PTV value
results_oar = {}
scale = np.arange(1.5, 0.4, -0.1)
pars = [6271, 2000, 5000, 6930, 3500]
ptv_vals = []
for ii in range(len(scale)):
    for jj in range(1, len(pars)):
        pars[jj] = scale[ii]*max_vals[jj]
    bayes2.set_pars(plan, pars)
    bayes2.calc_plan(plan, beam_set)
    objs = bayes2.get_objs(plan)
    results_oar[scale[ii]] = objs
    ptv = objs['PTV 4/7/20'][-1]['ResultValue']
    ptv_vals.append(ptv)
    print(f'Scale: {scale[ii]}, Max: {ptv}')
fpath = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\results\\'
np.save(fpath + 'oar_vals.npy', ptv_vals)
    
# # Look for percent difference with lowest PTV value
# results_ptv = {}
# scale = np.arange(1, 1.225, 0.025)
# pars = [6271, 2000, 5000, 6930, 3500]
# ptv_vals = []
# plt.figure()
# plt.xlabel('Scale Value')
# plt.ylabel('PTV Max')
# plt.title('Changing PTV Par')
# for ii in range(len(scale)):
#     pars[0] = scale[ii]*max_vals[0]
#     bayes2.set_pars(plan, pars)
#     bayes2.calc_plan(plan, beam_set)
#     objs = bayes2.get_objs(plan)
#     results_ptv[scale[ii]] = objs
#     ptv = objs['PTV 4/7/20'][-1]['ResultValue']
#     ptv_vals.append(ptv)
#     plt.plot(scale[ii], ptv, '.')
#     plt.draw()
# np.save(fpath + 'ptv_vals.npy', ptv_vals)
