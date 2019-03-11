#coding:utf-8

# A class of iir Low Shelving filter

# Check version
#  Python 3.6.4 on win32 (Windows 10)
#  numpy 1.14.0 
#  scipy 1.0.0
#  matplotlib  2.1.1


import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from HPF import *


class Class_IIR_LowShelving1(object):
    def __init__(self, fc=150, gain=42.0, slope=0.51, sampling_rate=48000):
        # design iir Low Shelving filter
        # initalize
        self.fc= fc # midpoint frequecny by unit is [Hz]
        self.sr= sampling_rate # sampling frequecny by unit is [Hz]
        self.gain= gain # amplification factor (magnification).   This must be > 0.0
        self.gain= np.sqrt(self.gain)
        self.slope= slope # shelf slope (S=1 for steepest slope)
        self.b, self.a= self.set_lowshelving()
        #print ('self.b,self.a', self.b, self.a)

    def filtering(self, xin):
    	# process of filtering
    	# input xin
    	# output filtered xin
        return signal.lfilter(self.b, self.a, xin)
        
    def f_show(self, worN=1024):
        # show frequency response
        wlist, fres = signal.freqz(self.b, self.a, worN=worN)
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        flist = wlist / ((2.0 * np.pi) / self.sr)
        plt.title('frequency response')
        ax1 = fig.add_subplot(111)
        
        plt.semilogx(flist, 20 * np.log10(abs(fres)), 'b')  # plt.plot(flist, 20 * np.log10(abs(fres)), 'b')
        
        plt.ylabel('Amplitude [dB]', color='b')
        plt.xlabel('Frequency [Hz]')
        
        ax2 = ax1.twinx()
        angles = np.unwrap(np.angle(fres))
        angles = angles / ((2.0 * np.pi) / 360.0)
        plt.semilogx(flist, angles, 'g')  # plt.plot(flist, angles, 'g')
        plt.ylabel('Angle(deg)', color='g')
        plt.grid()
        plt.axis('tight')
        plt.show()
        
    def f_show2(self, worN=1024):
        # compare low shelving with HPF
        # get HPF
        hpf=Class_HPF()
        wlist2, fres2 = signal.freqz(hpf.b, hpf.a, worN=worN)
        
        # show frequency response
        wlist, fres = signal.freqz(self.b, self.a, worN=worN)
        
        #
        fres3= fres2 * fres.T
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        flist = wlist / ((2.0 * np.pi) / self.sr)
        flist2 = wlist2 / ((2.0 * np.pi) / self.sr)
        plt.title('frequency response: blue lowshelving, green HPF, red combination')
        ax1 = fig.add_subplot(111)
        
        plt.semilogx(flist, 20 * np.log10(abs(fres)), 'b')  # plt.plot(flist, 20 * np.log10(abs(fres)), 'b')
        plt.semilogx(flist2, 20 * np.log10(abs(fres2)), 'g')  # plt.plot(flist, 20 * np.log10(abs(fres)), 'b')
        plt.semilogx(flist2, 20 * np.log10(abs(fres3)), 'r')  # plt.plot(flist, 20 * np.log10(abs(fres)), 'b')

        
        plt.ylabel('Amplitude [dB]', color='b')
        plt.xlabel('Frequency [Hz]')
        
        """
        ax2 = ax1.twinx()
        angles = np.unwrap(np.angle(fres))
        angles = angles / ((2.0 * np.pi) / 360.0)
        plt.semilogx(flist, angles, 'g')  # plt.plot(flist, angles, 'g')
        plt.ylabel('Angle(deg)', color='g')
        """
        plt.grid()
        plt.axis('tight')
        plt.show()

    def set_lowshelving(self,):
        
        omega= (self.fc / self.sr) * np.pi * 2.0
        sn= np.sin(omega)
        cs= np.cos(omega)
        alpha = sn / 2.0 * np.sqrt((self.gain + 1.0/self.gain) * (1.0/self.slope - 1.0) + 2.0)
        
        b=np.zeros(3) # umerator(bunsi)
        a=np.zeros(3) # denominator(bunbo)
        
        a[0]= ((self.gain + 1.0 ) + (self.gain - 1.0) * cs + 2.0 * np.sqrt(self.gain) * alpha)
        a[1]= -2.0 * (( self.gain - 1.0) + ( self.gain + 1.0) * cs)
        a[2]= ((self.gain + 1.0) + (self.gain - 1.0 ) * cs - 2.0 * np.sqrt( self.gain) * alpha )
        
        b[0]= self.gain * ((self.gain +1.0 ) - (self.gain - 1.0) * cs + 2.0 * np.sqrt(self.gain) * alpha )
        b[1]= 2.0 * self.gain * ((self.gain - 1.0) - (self.gain + 1.0 ) * cs)
        b[2]= self.gain * ((self.gain + 1.0) - (self.gain - 1.0) * cs - 2.0 * np.sqrt(self.gain) * alpha)
        
        b /= a[0]
        a /= a[0]
        
        return b, a


if __name__ == '__main__':
    
    # to get a low shelving filter as inverse filter against HPF that simulates radiation from mouth
    
    # low shelf sample 
    iir_LS1=Class_IIR_LowShelving1()  # fc=50, gain=410.0, slope=0.49)
    # draw frequency response
    iir_LS1.f_show2()
    
    
