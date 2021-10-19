# Import read toa from the write toa file under common io
# Directory is the my_shared_folder/EODP.../ed... ~ for a virtual machine based system

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
    check = False
    if counter > threesigma:                                            #Checking if three-sigma of the points are <0.01%
        check = True
    else:
        Check_pass = False
    return Check_pass


#1. Check for all bands that the differences with respect to the output TOA (ism_toa_isrf) are <0.01% for at least 3-sigma of the points.
#2. Check for all bands that the differences with respect to the output TOA (ism_toa_optical) are <0.01% for at least 3-sigma of the points.
band = [VNIR-0, VNIR-1, VNIR-2, VNIR-3]
nlines, ncolumns = 100, 150
# Read the toa's from the import - turn this into a loop
for i in range(len(band)):
    toa_input_1 = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_isrf' + str(band) + '.nc')
    toa_test_1 = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-ISM/test', 'ism_toa_isrf' + str(band) + '.nc')
    Check_passed_1 = threesigmacheck(toa_input_1, toa_test_1, nlines, ncolumns)
    toa_input_2 = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-L1B/input', 'ism_toa_optical' + str(band) + '.nc')
    toa_test_2 = readToa('/home/luss/my_shared_folder/EODP_TER_2021/EODP-TS-ISM/test', 'ism_toa_optical' + str(band) + '.nc')
    Check_passed_2 = threesigmacheck(toa_input_2, toa_test_2, nlines, ncolumns)
    if Check_passed_1 == True and Check_passed_2 == False:
        print("The differences with respect to the output are <0.01% for at least 3-sigma of the points.", str(band))

#3. What is the radiance to irradiance conversion factor for each band. What are the units of the TOA at this stage.
# See line 93 in opticalPhase.py at function Rad2Ird
# This is wrote in the rad2Irrad function -> should move to a seperate function

#4.a. Plot for all bands the System MTF across and along track (for the central pixels).
#4.b.  Report the MTF at the Nyquist frequency. Explain whether this is a decent or mediocre value and why.
# Fix plotMtf function in mtf.py; starting line 263
# Mtf values stored in a text document


#5. Explain the cause of the border effect introduced by the spatial filter (MTF)
# and what would be an appropriate solution (if any). How many pixel lines does it affect (roughly).

#6. Plot the TOA for all bands after the optical stage (with Panoply).
