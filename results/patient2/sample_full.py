"""Sample treatment plans for patient 2."""
import pandas as pd

import connect
import sample

# Inputs
repo_path = '\\\\client\\C$\\Users\Kelsey\\Dropbox (uwamath)\\autoray\\'
funcs_path = repo_path + 'results\\patient2\\funcs_full.csv'
goals_path = repo_path + 'results\\patient2\\goals.csv'
save_path = repo_path + 'results\\patient2\\'
max_iter = 1000
n_stop = 100

# Get RayStation objects
plan = connect.get_current('Plan')
beam_set = connect.get_current('BeamSet')

# Define functions and goals
funcs = sample.load_funcs(funcs_path)
goals = pd.read_csv(goals_path)
pars = sample.init_pars(funcs)
results = sample.init_results(goals)
stats = sample.init_stats()
roi_names = set(goals['Roi'])

# Sample treatment plans
n_success = 0
for ii in range(max_iter):
    print(f'Iteration: {ii}')
    if ii > 0:
        pars = sample.sample_pars(ii, funcs, pars)
        pars.to_pickle(save_path + 'pars_full.npy')
    print(f"Pars: {pars[pars['Sample']==ii]['DoseLevel'].values}")
    sample.set_pars(plan, pars)
    flag = sample.calc_plan(plan, beam_set, 'PTV 4/7/20', 6270, 99)
    if flag == 0:
        n_success += 1
        results = sample.get_results(plan, ii, flag, goals, results)
        stats = sample.get_stats(plan, ii, roi_names, stats)
        results.to_pickle(save_path + 'results_full.npy')
        stats.to_pickle(save_path + 'stats_full.npy')
    print(f'Success: {flag}, n_success: {n_success}\n')
    if n_success == n_stop:
        break
