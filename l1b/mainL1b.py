
# MAIN FUNCTION TO CALL THE L1B MODULE

from l1b.src.l1b import l1b

# Directory - this is the common directory for the execution of the E2E, all modules
auxdir = r'C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/eodp_students-master\/auxiliary'
#auxdir = '/home/luss/EODP/eodp/auxiliary'
indir = r'C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/EODP_TER_2021/EODP-TS-L1B/input'
outdir = r'C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/test_l1b'

# Initialise the ISM
myL1b = l1b(auxdir, indir, outdir)
myL1b.processModule()


