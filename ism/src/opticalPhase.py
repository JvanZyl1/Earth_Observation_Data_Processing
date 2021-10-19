from ism.src.initIsm import initIsm
from math import pi
from ism.src.mtf import mtf
from numpy.fft import fftshift, ifft2, fft2
import numpy as np
from common.io.writeToa import writeToa
from common.io.readIsrf import readIsrf
from scipy.interpolate import interp1d, interp2d
from common.plot.plotMat2D import plotMat2D
from common.plot.plotF import plotF
from scipy.signal import convolve2d
from common.src.auxFunc import getIndexBand
import matplotlib.pyplot as plt

class opticalPhase(initIsm):

    def __init__(self, auxdir, indir, outdir):
        super().__init__(auxdir, indir, outdir)

    def compute(self, sgm_toa, sgm_wv, band):
        """
        The optical phase is in charge of simulating the radiance
        to irradiance conversion, the spatial filter (PSF)
        and the spectral filter (ISRF).
        :return: TOA image in irradiances [mW/m2/nm],
                    with spatial and spectral filter
        """
        self.logger.info("EODP-ALG-ISM-1000: Optical stage")

        # Calculation and application of the ISRF
        # -------------------------------------------------------------------------------
        self.logger.info("EODP-ALG-ISM-1010: Spectral modelling. ISRF")
        toa = self.spectralIntegration(sgm_toa, sgm_wv, band)

        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [e-]")

        if self.ismConfig.save_after_isrf:
            saveas_str = self.globalConfig.ism_toa_isrf + band
            writeToa(self.outdir, saveas_str, toa)

        # Radiance to Irradiance conversion
        # -------------------------------------------------------------------------------
        self.logger.info("EODP-ALG-ISM-1020: Radiances to Irradiances")
        toa = self.rad2Irrad(toa,
                             self.ismConfig.D,
                             self.ismConfig.f,
                             self.ismConfig.Tr,
                             band)

        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [e-]")

        # Spatial filter
        # -------------------------------------------------------------------------------
        # Calculation and application of the system MTF
        self.logger.info("EODP-ALG-ISM-1030: Spatial modelling. PSF/MTF")
        myMtf = mtf(self.logger, self.outdir)
        Hsys = myMtf.system_mtf(toa.shape[0], toa.shape[1],
                                self.ismConfig.D, self.ismConfig.wv[getIndexBand(band)], self.ismConfig.f, self.ismConfig.pix_size,
                                self.ismConfig.kLF, self.ismConfig.wLF, self.ismConfig.kHF, self.ismConfig.wHF,
                                self.ismConfig.defocus, self.ismConfig.ksmear, self.ismConfig.kmotion,
                                self.outdir, band)
        #plt.show()

        # Apply system MTF
        toa = self.applySysMtf(toa, Hsys) # always calculated
        self.logger.debug("TOA [0,0] " +str(toa[0,0]) + " [e-]")
        # Write output TOA & plots
        # -------------------------------------------------------------------------------
        if self.ismConfig.save_optical_stage:
            saveas_str = self.globalConfig.ism_toa_optical + band

            writeToa(self.outdir, saveas_str, toa)

            title_str = 'TOA after the optical phase [mW/sr/m2]'
            xlabel_str='ACT'
            ylabel_str='ALT'
            plotMat2D(toa, title_str, xlabel_str, ylabel_str, self.outdir, saveas_str)

            idalt = int(toa.shape[0]/2)
            saveas_str = saveas_str + '_alt' + str(idalt)
            plotF([], toa[idalt,:], title_str, xlabel_str, ylabel_str, self.outdir, saveas_str)

        return toa

    def rad2Irrad(self, toa, D, f, Tr,band):
        """
        Radiance to Irradiance conversion
        :param toa: Input TOA image in radiances [mW/sr/m2]
        :param D: Pupil diameter [m]
        :param f: Focal length [m]
        :param Tr: Optical transmittance [-]
        :return: TOA image in irradiances [mW/m2]
        """
        omega = Tr*(pi/4)*((D/f)**2) #Equation page 34 of the reader : EODP-ALG-ISM-1020: This is the radiance to irradiance factor for each band
        toa = toa * omega #Radiance -> Irradiance

        #This is the testing procedure - to print the radiance to irradiance conversion factors
        PrintPlease = True
        listl = []
        if PrintPlease == True:                 #This is part of the test procedure,
            if str(band) == "VNIR-0":
                with open('/home/luss/my_shared_folder/EODP_Optical_Phase.txt', 'w') as f:
                    f.write("Start")
            with open('/home/luss/my_shared_folder/EODP_Optical_Phase.txt','r') as f:
                for line in f:
                    strip_lines=line.strip()
                    listli=strip_lines.split()
                    listl.append(listli)
            with open('/home/luss/my_shared_folder/EODP_Optical_Phase.txt', 'w') as f:
                f.write(str(listl))
                f.write("rad2Irrad_factor")
                f.write(str(band))
                f.write(str(omega))

        # TODO - DONE (8th October 2021)
        return toa

    def applySysMtf(self, toa, Hsys):
        """
        Application of the system MTF to the TOA
        :param toa: Input TOA image in irradiances [mW/m2]
        :param Hsys: System MTF
        :return: TOA image in irradiances [mW/m2]
        """
        # TODO
        GE = fft2(toa)              #1. Converting the TOA to the frequency domain.
        HE = fftshift(Hsys)         #2. Shift the zero frequency component to the center of the spectrum.
        FE = GE*HE                  #3. Multiply the shifted system MTF with the TOA in the frequency domain.
        toa_ft = ifft2(FE)          #4. The 2-dimensional inverse discrete Fourier Transform.
        #5. Check that the imaginary part is neglible and keep only the real part
        Imagi = np.imag(toa_ft)     #The imaginary components of the array
        neglible = 0.5              #Creates a neglible value upper bound
        infinitesmal = neglible * np.ones(np.shape(Imagi))   #Creates an array the size of "Imagi", with all components equal to "neglible"
        Imag_Bol = np.less_equal(Imagi, infinitesmal)
        Check = np.all(Imag_Bol)
        #if Check == True:
            #self.logger.info("The imaginary numbers of the TOA (after MTF is applied) are neglible, i.e. all are under", str(neglible))
        #else:
            #self.logger.info("The imaginary numbers of the TOA (after MTF is applied) aren't neglible, i.e. some are over", str(neglible))
        toa_ft = np.real(toa_ft)
        return toa_ft

    def spectralIntegration(self, sgm_toa, sgm_wv, band):
        """
        Integration with the ISRF to retrieve one band
        :param sgm_toa: Spectrally oversampled TOA cube 3D in irradiances [mW/m2]
        :param sgm_wv: wavelengths of the input TOA cube
        :param band: band
        :return: TOA image 2D in radiances [mW/m2]
        """
        # TODO - THIS NEEDS TO BE CHECKED
        # 1. Read the ISRF for it's band
        Isrf, wvisrf = readIsrf("/home/luss/my_shared_folder/eodp_students-master/auxiliary/isrf/","ISRF_" + str(band))
        # 2. Normalizing the ISRF
        Int_Isrf = np.sum(Isrf)
        Isrf_n = Isrf/Int_Isrf
        #####TO DO - align units of sgm_toa and sgm_wv such that they have the same units i.e. microns and nm
        #sigma_toa (100, 150, 600) i.e. (alt/nlines, act/ncolumns, wavelengths/nlambda)
        #nlines (alt, along-track) ~ 100, ncolumns (act, across-track) ~ 150, nlambda (spectral) ~ 600
        toa = np.zeros((sgm_toa.shape[0], sgm_toa.shape[1]))
        for ialt in range(sgm_toa.shape[0]):
            for iact in range(sgm_toa.shape[1]):
                cs = interp1d(sgm_wv,sgm_toa[ialt, iact, :], fill_value=(0, 0), bounds_error=False)
                toa_interp = cs(1000*wvisrf) #Check this is correct
                toa_in = toa_interp*Isrf_n
                toa[ialt, iact] = np.sum(toa_in)
        print("Shape of toa", np.shape(toa))
        return toa
