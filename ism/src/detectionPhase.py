from ism.src.initIsm import initIsm
import numpy as np
from common.io.writeToa import writeToa
from common.plot.plotMat2D import plotMat2D
from common.plot.plotF import plotF

import math
from scipy import stats

class detectionPhase(initIsm):

    def __init__(self, auxdir, indir, outdir):
        super().__init__(auxdir, indir, outdir)

        # Initialise the random see for the PRNU and DSNU
        np.random.seed(self.ismConfig.seed)


    def compute(self, toa, band):

        self.logger.info("EODP-ALG-ISM-2000: Detection stage")

        # Irradiance to photons conversion
        # -------------------------------------------------------------------------------
        self.logger.info("EODP-ALG-ISM-2010: Irradiances to Photons")
        area_pix = self.ismConfig.pix_size * self.ismConfig.pix_size # [m2]
        toa = self.irrad2Phot(toa, area_pix, self.ismConfig.t_int, self.ismConfig.wv[int(band[-1])])

        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [ph]")

        # Photon to electrons conversion
        # -------------------------------------------------------------------------------
        self.logger.info("EODP-ALG-ISM-2030: Photons to Electrons")
        toa = self.phot2Electr(toa, self.ismConfig.QE)

        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [e-]")

        if self.ismConfig.save_after_ph2e:
            saveas_str = self.globalConfig.ism_toa_e + band
            writeToa(self.outdir, saveas_str, toa)

        # PRNU
        # -------------------------------------------------------------------------------
        if self.ismConfig.apply_prnu:

            self.logger.info("EODP-ALG-ISM-2020: PRNU")
            toa = self.prnu(toa, self.ismConfig.kprnu)

            self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [e-]")

            if self.ismConfig.save_after_prnu:
                saveas_str = self.globalConfig.ism_toa_prnu + band
                writeToa(self.outdir, saveas_str, toa)

        # Dark-signal
        # -------------------------------------------------------------------------------
        if self.ismConfig.apply_dark_signal:

            self.logger.info("EODP-ALG-ISM-2020: Dark signal")
            toa = self.darkSignal(toa, self.ismConfig.kdsnu, self.ismConfig.T, self.ismConfig.Tref,
                                  self.ismConfig.ds_A_coeff, self.ismConfig.ds_B_coeff)

            self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [e-]")

            if self.ismConfig.save_after_ds:
                saveas_str = self.globalConfig.ism_toa_ds + band
                writeToa(self.outdir, saveas_str, toa)

        # Bad/dead pixels
        # -------------------------------------------------------------------------------
        if self.ismConfig.apply_bad_dead:

            self.logger.info("EODP-ALG-ISM-2050: Bad/dead pixels")
            toa = self.badDeadPixels(toa,
                               self.ismConfig.bad_pix,
                               self.ismConfig.dead_pix,
                               self.ismConfig.bad_pix_red,
                               self.ismConfig.dead_pix_red)


        # Write output TOA
        # -------------------------------------------------------------------------------
        if self.ismConfig.save_detection_stage:
            saveas_str = self.globalConfig.ism_toa_detection + band

            writeToa(self.outdir, saveas_str, toa)

            title_str = 'TOA after the detection phase [e-]'
            xlabel_str='ACT'
            ylabel_str='ALT'
            plotMat2D(toa, title_str, xlabel_str, ylabel_str, self.outdir, saveas_str)

            idalt = int(toa.shape[0]/2)
            saveas_str = saveas_str + '_alt' + str(idalt)
            plotF([], toa[idalt,:], title_str, xlabel_str, ylabel_str, self.outdir, saveas_str)
        return toa


    def irrad2Phot(self, toa, area_pix, tint, wv):
        """
        Conversion of the input Irradiances to Photons
        :param toa: input TOA in irradiances [mW/m2]
        :param area_pix: Pixel area [m2]
        :param tint: Integration time [s]
        :param wv: Central wavelength of the band [m]
        :return: Toa in photons
        """
        #TODO
        #toa = toa*(10**(-3)) #mW/m^2 -> W/m^2
        print("classical", toa[2,5])
        E_in = toa*area_pix*tint
        h = 6.62606896 * (10**(-34)) #Planck's constant
        C = 2.99792457 * (10**(8)) #Speed of light in a vacuum
        Ephotonk = (h*C)/wv
        toa_ph = E_in/Ephotonk
        print(toa_ph[2,5])
        return toa_ph

    def phot2Electr(self, toa, QE):
        """
        Conversion of photons to electrons
        :param toa: input TOA in photons [ph]
        :param QE: Quantum efficiency [e-/ph]
        :return: toa in electrons
        """
        #TODO
        toae = toa * QE
        #CHECK THE Ne < FWC
        print("toae", toae[2,5])
        return toae

    def badDeadPixels(self, toa,bad_pix,dead_pix,bad_pix_red,dead_pix_red):
        """
        Bad and dead pixels simulation
        :param toa: input toa in [e-]
        :param bad_pix: Percentage of bad pixels in the CCD [%]
        :param dead_pix: Percentage of dead pixels in the CCD [%]
        :param bad_pix_red: Reduction in the quantum efficiency for the bad pixels [-, over 1]
        :param dead_pix_red: Reduction in the quantum efficiency for the dead pixels [-, over 1]
        :return: toa in e- including bad & dead pixels
        """
        #TODO
        #1. Calculate the number of pixels affected
        n_col, n_row = np.shape(toa)
        size_toa = n_col*n_row
        n_pix_bad = n_col*(0.01*bad_pix)
        n_pix_dead = n_col*(0.01*dead_pix)
        #2. Assign these index locations with these relations, where toa_act is the number of
        #pixels in the across-track direction, and step bad and dead are the evenly distributed
        #steps (size_toa/n_pix).
        toa_act = n_row                                      #Number of pixels in the across-track direction
        step_bad = int(size_toa/n_pix_bad)
        step_dead = int(size_toa/n_pix_dead)
        idx_bad = range(5, toa_act, step_bad)
        idx_dead = range(0, toa_act, step_dead)

        #3. Apply the reduction factor to the DNS
        for i in range(len(idx_bad)):
            toa[:, idx_bad[i]] = toa[:, idx_bad[i]]*(1-bad_pix_red)   #CHECK THIS IS THE CORRECT WAY
        for j in range(len(idx_dead)):
            toa[:, idx_dead[j]] = toa[:, idx_dead[j]]*(1-dead_pix_red)  #CHECK THIS IS THE CORRECT WAY
        #4. Save to file (an ASCII txt file for example), the indexes, for validation purposes.

        index_bad = np.arange(5, toa_act, step_bad)
        index_dead = np.arange(0, toa_act, step_dead)
        index = list([index_bad, index_dead])
        PrintPlease = True
        if PrintPlease == True:                 #If you want the indexes set PrintPlease to False
            with open('/home/luss/my_shared_folder/EODP-ALG-ISM-2050.txt', 'w') as f:
                f.write("First index bad")
                f.write(str(index[0]))
                f.write("Second is index dead")
                f.write(str(index[1]))
        return toa

    def prnu(self, toa, kprnu):
        """
        Adding the PRNU effect
        :param toa: TOA pre-PRNU [e-]
        :param kprnu: multiplicative factor to the standard normal deviation for the PRNU
        :return: TOA after adding PRNU [e-]
        """
        #TODO
        n_col, n_row = np.shape(toa)
        n_act = n_col
        act = np.arange(1, n_act+1, 1)
        mean, sd = 0.0, 1.0
        f = np.zeros(np.shape(act))
        toa_a = np.zeros(np.shape(toa))
        for i in range(len(act)):
            prob_density = (1/(math.sqrt(2*np.pi)*sd)) * np.exp(-0.5*((act[i]-mean)/sd)**2)
            PRNU = prob_density*kprnu
            toa_a[i,:] = toa[i,:] * (1+PRNU)
        toa = toa_a
        return toa


    def darkSignal(self, toa, kdsnu, T, Tref, ds_A_coeff, ds_B_coeff):
        """
        Dark signal simulation
        :param toa: TOA in [e-]
        :param kdsnu: multiplicative factor to the standard normal deviation for the DSNU
        :param T: Temperature of the system
        :param Tref: Reference temperature of the system
        :param ds_A_coeff: Empirical parameter of the model 7.87 e-
        :param ds_B_coeff: Empirical parameter of the model 6040 K
        :return: TOA in [e-] with dark signal
        """
        #TODO

        n_col, n_row = np.shape(toa)
        n_act = n_col
        act = np.arange(1, n_act+1, 1)
        Sd = ds_A_coeff*((T/Tref)**3)*math.exp(-ds_B_coeff*(1/T - 1/Tref))                  #Constant component of the dark signal
        mean, sd = 0.0, 1.0
        f = np.zeros(np.shape(act))
        DS = np.zeros(np.shape(act))
        DSNU = np.zeros(np.shape(act))
        for i in range(len(act)):
            prob_density = (1/(math.sqrt(2*np.pi)*sd)) * np.exp(-0.5*((act[i]-mean)/sd)**2)
            d_s_n_u = abs(prob_density)*kdsnu
            DSNU[i] = d_s_n_u
            d_s = Sd*(1+d_s_n_u)                        #Total Dark Signal Changes Per Pixel
            DS[i] = d_s
            toa[i,:] = toa[i,:] + d_s    #Check this replaces
        return toa
