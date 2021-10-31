
# MAIN FUNCTION TO CALL THE L1C MODULE

from l1c.src.l1c import l1c

# Directory - this is the common directory for the execution of the E2E, all modules
auxdir = '/home/luss/my_shared_folder/eodp_students-master/auxiliary'
# GM dir + L1B dir
#indir = '/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1C/input/gm_alt100_act_150/,/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1C/input/l1b_output'
outdir = '/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1C/test'

#For E2E
indir = '/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-E2E/gm_out,/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-E2E/l1b_out'


# Initialise the ISM
myL1c = l1c(auxdir, indir, outdir)
myL1c.processModule()
