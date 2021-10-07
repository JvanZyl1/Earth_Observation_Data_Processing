# Import read toa from the write toa file under common io
# Directory is the my_shared_folder/EODP.../ed... ~ for a virtual machine based system

from common.io.writeToa import readToa
import numpy as np

# Read the toa's from the import - turn this into a loop
toa_0 = readToa('C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/test_l1b', 'l1b_toa_eq_VNIR-0.nc') #Perhaps .nc
toa_0_in = readToa(r'C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_VNIR-0.nc')
toa_1 = readToa('C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/test_l1b', 'l1b_toa_eq_VNIR-1.nc') #Perhaps .nc
toa_1_in = readToa(r'C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_VNIR-1.nc')
toa_2 = readToa('C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/test_l1b', 'l1b_toa_eq_VNIR-2.nc') #Perhaps .nc
toa_2_in = readToa(r'C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_VNIR-2.nc')
toa_3 = readToa('C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/test_l1b', 'l1b_toa_eq_VNIR-3.nc') #Perhaps .nc
toa_3_in = readToa(r'C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_VNIR-3.nc')
toa_all = [toa_0, toa_1, toa_2, toa_3]
toa_in = [toa_0_in, toa_1_in, toa_2_in, toa_3_in]

Check_pass = True #assumed it's passed
for i in range(len(toa_all)):                                           #Iterating over all the toa's bands
    N = 100*150                                                         #Calculating the number of pixels
    sum_ij = np.sum(toa_all[i])                                         #Finding the total toa
    mean = sum_ij/(N)                                                   #Calculating the mean of the values
    sigma = np.sqrt(np.sum(np.square(toa_all[i] - mean))/N**2)          #Calculating the standard deviation of the values (sigma)
    threesigma = 3*sigma                                                #Calculating three standard deviations
    difference = np.absolute(np.subtract(toa_all[i], toa_in[i]))        #Calculating the differences between the input and output
    percentage_diff_matrix = np.divide(difference, toa_0_in)*100        #Finding the percentage differences between the bands
    counter = 0
    for row in percentage_diff_matrix:
        for elem in row:
            if elem < 0.01:                                             #Checking if the differences are <0.01%
                counter += 1                                            #If so it counts it
    check = False
    if counter > threesigma:                                            #Checking if three-sigma of the points are <0.01%
        check = True
    else:
        Check_pass = False
if Check_pass == True:
    print("Question 1 is satisfied: for all the bands the differences with respect to the output are <0.01% for at least 3-sigma of the points.")


# Carry out the test tasks
#1. Find the Mean
    #mean = sum(radiance(i,j))/(no.pixels**2)
#2. Find the Standard Deivation
    #sigma = math.sqrt((sum(radiance(i,j) - mean))**2)/no.pixels
#3. Check sigma of each bank is <0.01% away for all points