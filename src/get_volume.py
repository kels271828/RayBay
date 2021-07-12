"""Get ROI names and volumes."""
import sys

repo_path = '\\\\client\\E$\\My Drive\\RayBay\\'
sys.path.append(repo_path + 'src\\')
import optimize

# Patient
patient_path = repo_path + 'results\\ZZ_MK_LULSBRT4544\\'

# Get ROI names and volumes
optimize.get_volumes(patient_path)
