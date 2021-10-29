from common.io.writeToa import readToa
import numpy as np

def threesigmacheck(toa_input, toa_test, nlines, ncolumns):
    N = nlines*ncolumns
    sum_ij = np.sum(toa_test)
    mean = sum_ij/N
    sigma = np.sqrt(np.sum(np.square(toa_test - mean))/N**2)
    threesigma = 3*sigma
    difference = np.absolute(np.subtract(toa_test, toa_input))
    percentage_diff_matrix = np.divide(difference, toa_input)*100
    counter = 0
    for row in percentage_diff_matrix:
        for elem in row:
            if elem < 0.01:                                             #Checking if the differences are <0.01%
                counter += 1                                            #If so it counts it
    Check_pass = False
    if counter > threesigma:                                            #Checking if three-sigma of the points are <0.01%
        Check_pass = True
    else:
        Check_pass = False
    return Check_pass

#1. Check for all bands that the differences with respect to the output TOA (ism_toa_isrf) are <0.01% for at least 3-sigma of the points.
#2. Check for all bands that the differences with respect to the output TOA (ism_toa_optical) are <0.01% for at least 3-sigma of the points.
band = ['VNIR-0', 'VNIR-1', 'VNIR-2', 'VNIR-3']
nlines, ncolumns = 100, 150
# Read the toa's from the import - turn this into a loop
for i in range(len(band)):
    toa_input_1 = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-ISM/output', 'ism_toa_' + band[i] + '.nc')
    toa_test_1 = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-ISM/test', 'ism_toa_' + band[i] + '.nc')
    Check_passed_1 = threesigmacheck(toa_input_1, toa_test_1, nlines, ncolumns)
    if Check_passed_1 == True:
        print("The differences with respect to the output are <0.01% for at least 3-sigma of the points.", band[i])
