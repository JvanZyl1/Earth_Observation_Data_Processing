# Import read toa from the write toa file under common io
# Directory is the my_shared_folder/EODP.../ed... ~ for a virtual machine based system

from common.io.writeToa import readToa
import numpy as np
import matplotlib.pyplot as plt

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
    percentage_diff_matrix = np.divide(difference, toa_in[i])*100        #Finding the percentage differences between the bands
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

isrf_0 = readToa('C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_isrf_VNIR-0.nc')
isrf_1 = readToa('C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_isrf_VNIR-1.nc')
isrf_2 = readToa('C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_isrf_VNIR-2.nc')
isrf_3 = readToa('C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_isrf_VNIR-3.nc')
isrf = [isrf_0, isrf_1, isrf_2, isrf_3]
isrf_row = []
for i in range(len(isrf)):
    b = np.sum(isrf[i], axis=0) #Sums columns of each array
    isrf_row.append(b)
print(np.shape(isrf_row)) #Checks that it's taking columns should be (4,150)
act = np.arange(1,151,1) #Creates a range of act pixels

l1b_toa_0 = readToa('C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1B/output', 'l1b_toa_VNIR-0.nc')
l1b_toa_1 = readToa('C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1B/output', 'l1b_toa_VNIR-1.nc')
l1b_toa_2 = readToa('C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1B/output', 'l1b_toa_VNIR-2.nc')
l1b_toa_3 = readToa('C:/Users/Jonathan van Zyl/Documents/BSC Aerospace Engineering/TU Delft_UC3M Yr.3/UC3M Exchange/Earth Observation and Data Processing/VM_shared_folder/EODP_TER_2021/EODP-TS-L1B/output', 'l1b_toa_VNIR-3.nc')
l1b_toa = [l1b_toa_0, l1b_toa_1, l1b_toa_2, l1b_toa_3]
l1b_row = []
for i in range(len(l1b_toa)):
    b = np.sum(l1b_toa[i], axis=0) #Sums columns of each array
    l1b_row.append(b)
print(np.shape(l1b_row)) #Checks that it's taking columns should be (4,150)


#Plotting of ism_toa_isrf and l1b_toa for each band with act_pixel vs. toa
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(act, isrf_row[0], 'tab:blue', act, l1b_row[0], 'tab:orange')
axs[0, 0].set_title('VNIR-0')
axs[0, 1].plot(act, isrf_row[1], 'tab:blue', act, l1b_row[1], 'tab:orange')
axs[0, 1].set_title('VNIR-1')
axs[1, 0].plot(act, isrf_row[2], 'tab:blue', act, l1b_row[2], 'tab:orange')
axs[1, 0].set_title('VNIR-2')
axs[1, 1].plot(act, isrf_row[3], 'tab:blue', act, l1b_row[3], 'tab:orange')
axs[1, 1].set_title('VNIR-3')

for ax in axs.flat:
    ax.set(xlabel='Act Pixel [-]', ylabel='TOA [mW/m2/sr')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

plt.show()
