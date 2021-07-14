"""Get ROI names and volumes."""
import sys

repo_path = '\\\\client\\E$\\My Drive\\RayBay\\'
sys.path.append(repo_path + 'src\\')
import optimize

# Patient
#patient_path = repo_path + 'results\\SBRT_lung_minsun\\'
#patient_path = repo_path + 'results\\ZZ_MK_LLungSBRT3778\\'
#patient_path = repo_path + 'results\\ZZ_MK_RLungSBRT4076\\'
#patient_path = repo_path + 'results\\ZZ_MK_RULungSBRT3796\\'
#patient_path = repo_path + 'results\\ZZ_MK_RLSBRT1931\\'
#patient_path = repo_path + 'results\\ZZ_MK_LULSBRT4544\\'
#patient_path = repo_path + 'results\\ZZ_MK_SBRTLL0924allviolated\\'
#patient_path = repo_path + 'results\\ZZ_MK_SBRTLL7289\\'
#patient_path = repo_path + 'results\\ZZ_MK_SBRTLLL8973\\'
#patient_path = repo_path + 'results\\ZZ_MK_SBRTRL7289\\'
patient_path = repo_path + 'results\\ZZ_MK_SBRTRUL_2928allviolate\\'

# Get ROI names and volumes
optimize.get_volumes(patient_path, 'roi_volumes.csv')
