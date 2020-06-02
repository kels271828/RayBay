"""Save results for specific ROIs."""
import utils

roi_names = ['SpinalCanal', 'Lungs', 'Lung_L', 'Lung_R', 'Heart',
             'Chestwall_L', 'Rib', 'PTV']
fpath = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\results\\'
utils.save_results(roi_names, fpath)
