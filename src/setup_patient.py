"""Get patient goals and functions.

Before running this script, run get_volumes.py to initialize goals.csv.
Next, fill in the columns `Type`, `GoalCriteria`, `DoseLevel (cGy)` and
`Volume (cm^3)`. Then calculate the column `Volume (%)` using the
equation 100*`Volume (cm^3)`/`RoiVolume (cm^3)`. Delete any rows not
included in the clinical goals.

"""
import os
import sys

repo_path = './'
sys.path.append(repo_path + 'src/')
import raybay

# Patient
patient_path = repo_path + 'results/ZZ_MK_LULSBRT4544/'

# Case
for case in ['approved', 'default', 'bayes']:
    if not os.path.exists(patient_path + case):
        os.makedirs(patient_path + case)
    raybay.create_goals(patient_path, case)
    raybay.create_funcs(patient_path, case)
