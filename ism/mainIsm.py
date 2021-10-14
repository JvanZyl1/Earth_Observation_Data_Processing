
# MAIN FUNCTION TO CALL THE ISM MODULE

from ism.src.ism import ism

# Directory - this is the common directory for the execution of the E2E, all modules
#auxdir = r'C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/eodp_students-master/auxiliary'
#indir = 'C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-ISM/input/gradient_alt100_act150'
#outdir = 'C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-ISM/test/'

#auxdir = '/home/luss/EODP/eodp/auxiliary'
auxdir = '/home/luss/my_shared_folder/eodp_students-master/auxiliary'
#indir = '/home/luss/my_shared_folder/sgm_out/gradient_alt100_act150/' # small scene
indir = '/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-ISM/input/gradient_alt100_act150/'
#outdir = '/home/luss/my_shared_folder/ism_out/'
outdir = '/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-ISM/test/'

# Initialise the ISM
myIsm = ism(auxdir, indir, outdir)
myIsm.processModule()
