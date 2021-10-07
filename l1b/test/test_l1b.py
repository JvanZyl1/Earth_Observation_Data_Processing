import nc as nc
import matplotlib.pyplot as plt
import numpy as np
import os

# Import read toa from the write toa file under common io
# Directory is the my_shared_folder/EODP.../ed...

from common.io.writeToa import readToa
import numpy as np
import math

# Read the toa's from the import - turn this into a loop
toa_0 = readToa('/home/luss/my_shared_folder/test_l1b', 'l1b_toa_eq_VNIR-0.nc') #Perhaps .nc
toa_0_in = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_VNIR-0.nc')
toa_1 = readToa('/home/luss/my_shared_folder/test_l1b', 'l1b_toa_eq_VNIR-1.nc')
toa_1_in = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_VNIR-1.nc')
toa_2 = readToa('/home/luss/my_shared_folder/test_l1b', 'l1b_toa_eq_VNIR-2.nc')
toa_2_in = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_VNIR-2.nc')
toa_3 = readToa('/home/luss/my_shared_folder/test_l1b', 'l1b_toa_eq_VNIR-3.nc')
toa_3_in = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_VNIR-3.nc')
toa_all = [toa_0, toa_1, toa_2, toa_3]
toa_in = [toa_0_in, toa_1_in, toa_2_in, toa_3_in]

Check_pass = True #assumed it's passed
for i in range(len(toa_all)):
    N = 100*150
    sum_ij = np.sum(toa_all[i])
    mean = sum_ij/(N)
    sigma = np.sqrt(np.sum(np.square(toa_all[i] - mean))/N**2)
    threesigma = 3*sigma
    difference = np.absolute(np.subtract(toa_all[i], toa_in[i]))
    percentage_diff_matrix = np.divide(difference, toa_0_in)*100
    counter = 0
    for row in percentage_diff_matrix:
        for elem in row:
            if elem < 0.01:
                counter += 1
    check = False
    if counter > threesigma:
        check = True
    else:
        Check_pass = False
if Check_pass == True:
    print("Question 1 is satisfied")

# Carry out the test tasks
#1. Find the Mean
    #mean = sum(radiance(i,j))/(no.pixels**2)
#2. Find the Standard Deivation
    #sigma = math.sqrt((sum(radiance(i,j) - mean))**2)/no.pixels
#3. Check sigma of each bank is <0.01% away for all points
#Link this python file to the main file so that it will run properly.
