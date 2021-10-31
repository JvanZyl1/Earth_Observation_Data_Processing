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


plots = False
if plots == False:
    exit(0)
'''PART 2: Plotting'''
act = np.arange(1,151,1) #Creates a range of act pixels
band = ['VNIR-0', 'VNIR-1', 'VNIR-2', 'VNIR-3']
nlines, ncolumns = 100, 150
isrf_row, l1b_row = [], []
# Read the toa's from the import - turn this into a loop
for i in range(len(band)):
    toa_isrf = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-ISM/output', 'ism_toa_' + band[i] + '.nc')
    #toa_restored = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1B/output', 'l1b_toa_' + band[i] + '.nc')
    toa_restored = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-ISM/test', 'ism_toa_' + band[i] + '.nc')
    isrf_row.append(np.sum(toa_isrf, axis=0))
    l1b_row.append(np.sum(toa_restored, axis=0))


Q3 = True #If you want the Q3
if Q3 == True:
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot( act, l1b_row[0], 'tab:orange')
    axs[0, 0].set_title('VNIR-0')
    axs[0, 1].plot( act, l1b_row[1], 'tab:orange')
    axs[0, 1].set_title('VNIR-1')
    axs[1, 0].plot( act, l1b_row[2], 'tab:orange')
    axs[1, 0].set_title('VNIR-2')
    axs[1, 1].plot( act, l1b_row[3], 'tab:orange')
    axs[1, 1].set_title('VNIR-3')

    for ax in axs.flat:
        ax.set(xlabel='Act Pixel [-]', ylabel='TOA [mW/m2/sr]')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

    plt.show()



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
    ax.set(xlabel='Act Pixel [-]', ylabel='TOA [mW/m2/sr]')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

plt.show()
