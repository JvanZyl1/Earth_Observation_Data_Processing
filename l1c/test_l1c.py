from common.io.writeToa import readToa
import numpy as np
import matplotlib.pyplot as plt

def threesigmacheck(toa_input, toa_test, nlines, ncolumns):
    N = nlines*ncolumns
    sum_ij = np.sum(toa_test)
    mean = sum_ij/N
    sigma = np.sqrt(np.sum(np.square(toa_test - mean))/N**2)
    threesigma = 3*sigma
    difference = np.absolute(np.subtract(toa_test, toa_input))
    #To avoid zero division error - this process takes that if the toa is zero then it takes the value difference value
    percentage_diff_matrix = np.zeros(np.shape(difference))
    for ii in range(len(difference)):
        if toa_input[ii] != 0:
            percentage_diff_matrix[ii] = np.divide(difference[ii], toa_input[ii])*100
        else:
            percentage_diff_matrix[ii] = 0
    counter = 0
    for elem in percentage_diff_matrix:
        if elem < 0.01:                                             #Checking if the differences are <0.01%
            counter += 1                                            #If so it counts it
    Check_pass = False
    if counter > threesigma:                                            #Checking if three-sigma of the points are <0.01%
        Check_pass = True
    else:
        Check_pass = False
    return Check_pass


band = ['VNIR-0', 'VNIR-1', 'VNIR-2', 'VNIR-3']
nlines, ncolumns = 100, 150
# Read the toa's from the import - turn this into a loop
for i in range(len(band)):
    toa_input = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1C/output', 'l1c_toa_' + band[i] + '.nc')
    toa_test = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1C/test', 'l1c_toa_' + band[i] + '.nc')
    # The output reference L1C TOA has some points with negative radiances which are the result of an extrapolation. If you want, you can set these negative radiances to zero.
    toa_input_a = np.where(toa_input > 0, toa_input, 0)
    # NOTE: The ordering of the TOAs is not the same in all executions, therefore, to check this requirement you need to sort the TOAs in your results and in the reference data.
    toa_input_1 = np.sort(toa_input_a)
    toa_test_1 = np.sort(toa_test)
    Check_passed = threesigmacheck(toa_input_1, toa_test_1, nlines, ncolumns)
    if Check_passed == True:
        print("The differences with respect to the output are <0.01% for at least 3-sigma of the points.", band[i])
    else:
        exit(1)
