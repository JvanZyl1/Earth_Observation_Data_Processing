
from ism.src.initIsm import initIsm
import numpy as np
from common.plot.plotMat2D import plotMat2D
from common.plot.plotF import plotF

class videoChainPhase(initIsm):

    def __init__(self, auxdir, indir, outdir):
        super().__init__(auxdir, indir, outdir)

    def compute(self, toa, band):
        self.logger.info("EODP-ALG-ISM-3000: Video Chain")

        # Electrons to Voltage - read-out & amplification
        # -------------------------------------------------------------------------------
        self.logger.info("EODP-ALG-ISM-3010: Electrons to Voltage – Read-out and Amplification")
        toa = self.electr2Volt(toa,
                         self.ismConfig.OCF,
                         self.ismConfig.ADC_gain)

        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [V]")

        # Digitisation
        # -------------------------------------------------------------------------------
        self.logger.info("EODP-ALG-ISM-3020: Voltage to Digital Numbers – Digitisation")
        toa = self.digitisation(toa,
                          self.ismConfig.bit_depth,
                          self.ismConfig.min_voltage,
                          self.ismConfig.max_voltage)

        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [DN]")

        # Plot
        if self.ismConfig.save_vcu_stage:
            saveas_str = self.globalConfig.ism_toa_vcu + band
            title_str = 'TOA after the VCU phase [DN]'
            xlabel_str='ACT'
            ylabel_str='ALT'
            plotMat2D(toa, title_str, xlabel_str, ylabel_str, self.outdir, saveas_str)

            idalt = int(toa.shape[0]/2)
            saveas_str = saveas_str + '_alt' + str(idalt)
            plotF([], toa[idalt,:], title_str, xlabel_str, ylabel_str, self.outdir, saveas_str)

        return toa

    def electr2Volt(self, toa, OCF, gain_adc):
        """
        Electron to Volts conversion.
        Simulates the read-out and the amplification
        (multiplication times the gain).
        :param toa: input toa in [e-]
        :param OCF: Output Conversion factor [V/e-]
        :param gain_adc: Gain of the Analog-to-digital conversion [-]
        :return: output toa in [V]
        """
        #TODO
        toa = toa*OCF*gain_adc
        PrintPlease = True
        convd = OCF*gain_adc
        if PrintPlease == True:                 #If you want the indexes set PrintPlease to False
            with open('/home/luss/my_shared_folder/EODP-CONVERSION_VID.txt', 'a') as f:
                a = ("\n", "Electron to Volts:", str(convd))
                str1 = ''
                b = str1.join(a)
                f.write(str(b))
                f.close()

        return toa

    def digitisation(self, toa, bit_depth, min_voltage, max_voltage):
        """
        Digitisation - conversion from Volts to Digital counts
        :param toa: input toa in [V]
        :param bit_depth: bit depth
        :param min_voltage: minimum voltage
        :param max_voltage: maximum voltage
        :return: toa in digital counts
        """
        #TODO
        toa_dn = np.round((toa/(max_voltage - min_voltage))*((2**bit_depth)-1))
        #https://www.google.com/search?q=bit+depth+formula&rlz=1C1CHBF_enNL923NL923&sxsrf=AOaemvIG6J57vAbOWelIekv9--G5X7FFAw:1635517427491&tbm=isch&source=iu&ictx=1&fir=bqDCUZcFtL2OSM%252CFjH_GSO4ypsBCM%252C_%253BPS1PCKuUPa7lHM%252C18SK0X0PLV8MRM%252C_%253BKYWgr5LLp-m9CM%252C7Jv1i5P3RL3PWM%252C_%253BNXUlAr2UEcS9JM%252C1sPuVNfinCtyXM%252C_%253B5tQC_qJhnRbaXM%252CmtafdBKCsWqnSM%252C_%253BtK5Qy7bz3KCPNM%252CmtafdBKCsWqnSM%252C_&vet=1&usg=AI4_-kRevuYinqvOrR266yxTPmgVDZPKYg&sa=X&ved=2ahUKEwiS5Kea6e_zAhULNhoKHRF8CAoQ9QF6BAgyEAE&biw=1536&bih=722&dpr=1.25#imgrc=PS1PCKuUPa7lHM
        bit_depth_DN = 2**(bit_depth) - 1
        counter = 0
        for i in range(np.shape(toa_dn)[0]):
            for j in range(np.shape(toa_dn)[1]):
                if toa_dn[i,j] > bit_depth_DN:
                    toa_dn[i,j] = bit_depth_DN
                    counter += 1
                if toa_dn[i,j] < 0:
                    toa_dn[i,j] = 0
        toa = toa_dn
        perc = counter/(np.size(toa_dn))

        PrintPlease2 = True
        if PrintPlease2 == True:                 #If you want the indexes set PrintPlease to False
            with open('/home/luss/my_shared_folder/EODP-Saturated.txt', 'a') as f:
                a = ("\n", "Percentage of Saturated pixels", str(perc))
                str1 = ''
                b = str1.join(a)
                f.write(str(b))
                f.close()


        convf = np.round((1/(max_voltage - min_voltage))*((2**bit_depth)-1))
        PrintPlease = True
        if PrintPlease == True:                 #If you want the indexes set PrintPlease to False
            with open('/home/luss/my_shared_folder/EODP-CONVERSION_VID.txt', 'a') as f:
                a = ("\n", "Volts to DN:", str(convf))
                str1 = ''
                b = str1.join(a)
                f.write(str(b))
                f.close()

        return toa

