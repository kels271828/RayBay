"""Save results for specific ROIs."""
import raystation

roi_names = ['SpinalCanal', 'Lungs', 'Lung_L', 'Lung_R', 'Heart',
             'Chestwall_L', 'Rib', 'PTV']
fpath = '\\\\client\\C$\\Users\\Kelsey\\Dropbox (uwamath)\\autoray\\results\\'
raystation.save_results(roi_names, fpath)
