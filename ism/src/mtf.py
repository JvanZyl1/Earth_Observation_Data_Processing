from math import pi
from config.ismConfig import ismConfig
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import j1
from numpy.matlib import repmat
from common.io.readMat import writeMat
from common.plot.plotMat2D import plotMat2D
from scipy.interpolate import interp2d
from numpy.fft import fftshift, ifft2
import os

class mtf:
    """
    Class MTF. Collects the analytical modelling of the different contributions
    for the system MTF
    """
    def __init__(self, logger, outdir):
        self.ismConfig = ismConfig()
        self.logger = logger
        self.outdir = outdir

    def system_mtf(self, nlines, ncolumns, D, lambd, focal, pix_size,
                   kLF, wLF, kHF, wHF, defocus, ksmear, kmotion, directory, band):
        """
        System MTF
        :param nlines: Lines of the TOA
        :param ncolumns: Columns of the TOA
        :param D: Telescope diameter [m]
        :param lambd: central wavelength of the band [m]
        :param focal: focal length [m]
        :param pix_size: pixel size in meters [m]
        :param kLF: Empirical coefficient for the aberrations MTF for low-frequency wavefront errors [-]
        :param wLF: RMS of low-frequency wavefront errors [m]
        :param kHF: Empirical coefficient for the aberrations MTF for high-frequency wavefront errors [-]
        :param wHF: RMS of high-frequency wavefront errors [m]
        :param defocus: Defocus coefficient (defocus/(f/N)). 0-2 low defocusing
        :param ksmear: Amplitude of low-frequency component for the motion smear MTF in ALT [pixels]
        :param kmotion: Amplitude of high-frequency component for the motion smear MTF in ALT and ACT
        :param directory: output directory
        :return: mtf
        """

        self.logger.info("Calculation of the System MTF")

        # Calculate the 2D relative frequencies
        self.logger.debug("Calculation of 2D relative frequencies")
        fn2D, fr2D, fnAct, fnAlt, w_inv = self.freq2d(nlines, ncolumns, D, lambd, focal, pix_size)

        # Diffraction MTF
        self.logger.debug("Calculation of the diffraction MTF")
        Hdiff = self.mtfDiffract(fr2D)

        # Defocus
        Hdefoc = self.mtfDefocus(fr2D, defocus, focal, D)

        # WFE Aberrations
        Hwfe = self.mtfWfeAberrations(fr2D, lambd, kLF, wLF, kHF, wHF)

        # Detector
        Hdet  = self. mtfDetector(fn2D)

        # Smearing MTF
        Hsmear = self.mtfSmearing(fnAlt, ncolumns, ksmear)

        # Motion blur MTF
        Hmotion = self.mtfMotion(fn2D, kmotion)

        # Calculate the System MTF
        self.logger.debug("Calculation of the Sysmtem MTF by multiplying the different contributors")
        Hsys = Hdiff*Hwfe*Hdefoc*Hdet*Hsmear*Hmotion                                           #TODO

        #Calculate the Nyquist frequency ---- SELF FUNCTION MADE
        fNyq = self.fNyuq(pix_size)
        self.logger.info("Calculation of the Nyquist frequency")

        # Plot cuts ACT/ALT of the MTF
        self.plotMtf(Hdiff, Hdefoc, Hwfe, Hdet, Hsmear, Hmotion, Hsys, nlines, ncolumns, fnAct, fnAlt, directory, band, fNyq, w_inv)


        return Hsys

    def freq2d(self,nlines, ncolumns, D, lambd, focal, w):
        """
        Calculate the relative frequencies 2D (for the diffraction MTF)
        :param nlines: Lines of the TOA
        :param ncolumns: Columns of the TOA
        :param D: Telescope diameter [m]
        :param lambd: central wavelength of the band [m]
        :param focal: focal length [m]
        :param w: pixel size in meters [m]
        :return fn2D: normalised frequencies 2D (f/(1/w))
        :return fr2D: relative frequencies 2D (f/(1/fc))
        :return fnAct: 1D normalised frequencies 2D ACT (f/(1/w))
        :return fnAlt: 1D normalised frequencies 2D ALT (f/(1/w))
        """
        print(D, "D", lambd, "lamd", focal, "focal")
        infinitesimal = 10**(-10)
        fstepAlt = 1/nlines/w
        fstepAct = 1/ncolumns/w
        fAlt = np.arange(-1/(2*w), 1/(2*w)-infinitesimal, fstepAlt)
        fAct = np.arange(-1/(2*w), 1/(2*w)-infinitesimal, fstepAct) #Reader page 63
        eps_cutoff = D/(lambd*focal) # Optical cutoff frequency
        fr_factor = eps_cutoff#(1/w)/eps_cutoff # As said in reader, shown below:
        #For the 2D relative frequencies, multiply by a factor (1/w)/ξc. This way we get the relative frequencies but with the same sampling than the fn2D calculated above.
        w_inv = 1/w
        fnAct = np.divide(fAct, w_inv) # The normalized frequencies 1D ACT
        fnAlt = np.divide(fAlt, w_inv) # The normalized frequencies 1D ALT
        frAct = np.divide(fAct, fr_factor) #The relative frequencies 1D ACT
        frAlt = np.divide(fAlt, fr_factor) #The relative frequencies 1D ALT

        [fnAltxx, fnActxx] = np.meshgrid(fnAlt, fnAct, indexing='ij')  # Please use ‘ij’ indexing or you will get the transpose
        fn2D = np.sqrt(fnAltxx * fnAltxx + fnActxx * fnActxx)

        [frAltxx, frActxx] = np.meshgrid(frAlt, frAct, indexing='ij')
        fr2D = np.sqrt(frAltxx * frAltxx + frActxx * frActxx)
        #TODO - DONE (8th October 2021)
        return fn2D, fr2D, fnAct, fnAlt, w_inv

    def mtfDiffract(self,fr2D):
        """
        Optics Diffraction MTF
        :param fr2D: 2D relative frequencies (f/fc), where fc is the optics cut-off frequency
        :return: diffraction MTF
        """
        #TODO
        #Hdiff = 1 - fr2D
        Hdiffb = []
        #print("1")
        #print(fr2D)
        #print("2")
        for i in range(len(fr2D)):
            Hdiff_a = 2/pi * (np.arccos(fr2D[i]) - fr2D[i]*(1 - np.power(fr2D[i],2))**0.5)
            Hdiffb.append(Hdiff_a)
        Hdiff = np.array(Hdiffb)
        #print(np.shape(Hdiff))

        #if Hdiff[fr2D*fr2D>1] = 0:                                 #This needs to be correctly implemented
        #    print("MTF check is completed")
        return Hdiff


    def mtfDefocus(self, fr2D, defocus, focal, D):
        """
        Defocus MTF
        :param fr2D: 2D relative frequencies (f/fc), where fc is the optics cut-off frequency
        :param defocus: Defocus coefficient (defocus/(f/N)). 0-2 low defocusing
        :param focal: focal length [m]
        :param D: Telescope diameter [m]
        :return: Defocus MTF
        """
        coe = pi*defocus
        minus = np.add(fr2D, -1)
        b = np.multiply(fr2D, minus)
        x = np.multiply(coe, b)
        #x = np.multiply((pi*defocus), (np.multiply(fr2D, (np.add(fr2D,-1))))) # This is for a more consisess equation
        #Bes_J1 = (np.divide(x,2)) - (np.divide((np.power(x,3)),16)) + (np.divide((np.power(x,5)),384)) - (np.divide((np.power(x,7)),18432))
        #Note that f/D hasn't been used
        Bes_J1 = j1(x)
        Hdefoc = 2*np.divide(Bes_J1, x) #Pages 51-52
        #print(np.shape(Hdefoc))
        return Hdefoc

    def mtfWfeAberrations(self, fr2D, lambd, kLF, wLF, kHF, wHF):
        """
        Wavefront Error Aberrations MTF
        :param fr2D: 2D relative frequencies (f/fc), where fc is the optics cut-off frequency
        :param lambd: central wavelength of the band [m]
        :param kLF: Empirical coefficient for the aberrations MTF for low-frequency wavefront errors [-]
        :param wLF: RMS of low-frequency wavefront errors [m]
        :param kHF: Empirical coefficient for the aberrations MTF for high-frequency wavefront errors [-]
        :param wHF: RMS of high-frequency wavefront errors [m]
        :return: WFE Aberrations MTF
        """
        inp = -fr2D*(1-fr2D)*((kLF*((wLF/lambd)**2)) + (kHF*((wHF/lambd)**2)))
        Hwfe = np.exp(inp)
        #print(np.shape(Hwfe))
        #TODO
        return Hwfe

    def mtfDetector(self,fn2D):
        """
        Detector MTF
        :param fnD: 2D normalised frequencies (f/(1/w))), where w is the pixel width
        :return: detector MTF
        """
        #TODO
        Hdet = np.abs(np.sinc(fn2D))
        #print(np.shape(Hdet))
        return Hdet

    def mtfSmearing(self, fnAlt, ncolumns, ksmear):
        """
        Smearing MTF
        :param ncolumns: Size of the image ACT
        :param fnAlt: 1D normalised frequencies 2D ALT (f/(1/w))
        :param ksmear: Amplitude of low-frequency component for the motion smear MTF in ALT [pixels]
        :return: Smearing MTF
        """
        # Calculate the 1D MTF in the ALT direction (using fn2D, size n-lines), and then repeat it in the ACT direction with a repmat
        Hsmear_alt = np.sinc(ksmear*fnAlt)
        Hsmeart = repmat(Hsmear_alt, ncolumns, 1)                    #Note for self, check repmat is working, as first time using function.
        Hsmear = np.transpose(Hsmeart)  #THIS HAS BEEN CHANGED TO ALTER ARRAY SHAPE - CHECK !!!!
        #print(np.shape(Hsmear))
        #TODO
        return Hsmear

    def mtfMotion(self, fn2D, kmotion):
        """
        Motion blur MTF
        :param fnD: 2D normalised frequencies (f/(1/w))), where w is the pixel width
        :param kmotion: Amplitude of high-frequency component for the motion smear MTF in ALT and ACT
        :return: detector MTF
        """
        unp = np.multiply(kmotion, fn2D)
        Hmotion = np.sinc(unp)
        #print(np.shape(Hmotion))
        #TODO
        return Hmotion

    def fNyuq(self, w):  #Self-Made frequency
        """
        Nyquist frequency calculation
        :param w: pixel size [m]
        :return: Nyquist freqency [-]
        """
        fNyq = 1/(2*w)
        return fNyq

    def plotMtf(self,Hdiff, Hdefoc, Hwfe, Hdet, Hsmear, Hmotion, Hsys, nlines, ncolumns, fnAct, fnAlt, directory, band, fNyq, w_inv):
        """
        Plotting the system MTF and all of its contributors
        :param Hdiff: Diffraction MTF
        :param Hdefoc: Defocusing MTF
        :param Hwfe: Wavefront electronics MTF
        :param Hdet: Detector MTF
        :param Hsmear: Smearing MTF
        :param Hmotion: Motion blur MTF
        :param Hsys: System MTF
        :param nlines: Number of lines in the TOA
        :param ncolumns: Number of columns in the TOA
        :param fnAct: normalised frequencies in the ACT direction (f/(1/w))
        :param fnAlt: normalised frequencies in the ALT direction (f/(1/w))
        :param directory: output directory
        :param band: band
        :return: N/A
        """
        H = np.array([Hdiff, Hdefoc, Hwfe, Hdet, Hsmear, Hmotion, Hsys])
        #print(np.shape(H)): 7 H values of shape (2, 801) == (nlines, ncolumns)
        #System MTF - Slice ACT & ACT
        # (Alt, Act) - Alt is rows, Act is columns
        mid_row = round(nlines / 2)
        mid_col = round(ncolumns /2)
        nyqnorm = fNyq/w_inv
        fna_mid = fnAlt
        H_mid_act = []
        H_mid_alt = []
        for i in range(len(H)):
            Hi = H[i]
            H_mid_alt.append(Hi[mid_row,:])
            H_mid_act.append(Hi[:,mid_col])

        #Plot of the ACT Slice
        plt.figure(3)
        plt.plot(fna_mid, H_mid_act[0], fna_mid, H_mid_act[1], fna_mid, H_mid_act[2],fna_mid, H_mid_act[1], fna_mid,
                 H_mid_act[3], fna_mid, H_mid_act[4], fna_mid, H_mid_act[5], fna_mid, H_mid_act[6])
        plt.xlim(0, 0.55)
        plt.axvline(x = nyqnorm, color = 'b', label = 'axvline - full height')
        #plt.vlines(fNyq, 0, 1, colors = 'black', linestyles = 'dashed')
        plt.legend([ "Nyquist", "Hdiff", "Hdefoc", "Hwfe", "Hdet", "Hsmear", "Hmotion", "Hsys"])
        plt.xlabel("Spatial frequencies f/(1/w) [-]")
        plt.ylabel("MTF")
        c = ("System MTF - slice ACT", band)
        plt.title(c)
        saveas_str = ("ism_act_slice", band)
        plt.savefig(directory + os.path.sep + saveas_str[0] + saveas_str[1] +'.png')
        plt.close()

        #Plot of the ALT Slice
        plt.figure(1)
        plt.plot(fnAct, H_mid_alt[0], fnAct, H_mid_alt[1], fnAct, H_mid_alt[2], fnAct, H_mid_alt[3], fnAct,
                 H_mid_alt[4], fnAct, H_mid_alt[5], fnAct, H_mid_alt[6])
        plt.xlim(0, 0.55)
        plt.axvline(x = nyqnorm, color = 'b', label = 'axvline - full height')
        plt.legend([ "Nyquist", "Hdiff", "Hdefoc", "Hwfe", "Hdet", "Hsmear", "Hmotion", "Hsys"])
        plt.xlabel("Spatial frequencies f/(1/w) [-]")
        plt.ylabel("MTF")
        b = ("System MTF - slice ALT", band)
        plt.title(b)
        saveas_str = ("ism_alt_slice", band)
        plt.savefig(directory + os.path.sep + saveas_str[0] + saveas_str[1] +'.png')
        plt.close()

        #System MTF - Contour
        plt.figure(3)
        X,Y = np.meshgrid(fnAct, fnAlt)
        Z = Hsys
        plt.contourf(X, Y, Z)
        plt.colorbar()
        plt.xlabel("ACT")
        plt.ylabel("ALT")
        a = ("System MTF for", band)
        plt.title(str(a))
        saveas_str = ("ism_system_mtf_contour", band)
        plt.savefig(directory + os.path.sep + saveas_str[0] + saveas_str[1] +'.png')
        plt.close(3)

        #Saving the MTF values at the Nyquist Frequency
