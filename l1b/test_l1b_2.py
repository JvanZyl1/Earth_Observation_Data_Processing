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
    for i in range(0,nlines, 1):
        for j in range(0, ncolumns, 1):
            percentage_diff_matrix[i,j] = np.divide(difference[i,j], toa_input[i,j])*100
    counter = 0
    for row in percentage_diff_matrix:
        for elem in row:
            if elem < 0.01:                                            #Checking if the differences are <0.01%
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
    toa_input_1 = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1B/output', 'l1b_toa_eq_' + band[i] + '.nc')
    toa_test_1 = readToa('/home/luss/my_shared_folder/test_l1b', 'l1b_toa_eq_' + band[i] + '.nc')
    Check_passed_1 = threesigmacheck(toa_input_1, toa_test_1, nlines, ncolumns)
    if Check_passed_1 == True:
        print("The differences with respect to the output are <0.01% for at least 3-sigma of the points.", str(band[i]))
    else:
        print("Differences are too great")
