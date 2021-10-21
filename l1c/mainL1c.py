
# MAIN FUNCTION TO CALL THE L1C MODULE

from l1c.src.l1c import l1c

# Directory - this is the common directory for the execution of the E2E, all modules
#auxdir = '/home/luss/EODP/eodp/auxiliary'
auxdir = '/home/luss/my_shared_folder/eodp_students-master/auxiliary'
#auxdir = r'C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/eodp_students-master\/auxiliary'
# GM dir + L1B dir
indir = '/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1C/input/gm_alt100_act_150/,/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1C/input/l1b_output'
#indir = '/home/luss/my_shared_folder/gm_out/gm_alt100_act_150/,/home/luss/EODP/eodp/l1b/test/ut02/output'
outdir = '/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1C/test'
#outdir = 'C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1C/output'


# Initialise the ISM
myL1c = l1c(auxdir, indir, outdir)
myL1c.processModule()
