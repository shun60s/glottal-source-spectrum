#coding:utf-8

#  return waveform of which phase is as same as specified ref_phase waveform, keeping original amplitude.
#  by FFT and inverse FFT
#
#　与えられた区間のデータをFFT分析により　周波数と位相特性を求める。
#  周波数特性はそのままで、ref_phaseの波形の位相特性と同じにしたものを返す。

import numpy as np
from scipy import signal
from scipy import fftpack
from matplotlib import pyplot as plt

# Check version
#  Python 3.6.4, 64bit on Win32 (Windows 10)
#  numpy (1.14.0)
#  scipy (1.0.0)

def fft_ana(xin, ref_phase, sr, comment0='', show=False):
    # check
    if len(xin) < 2:
        print ('warning: xin is too small to compute.')
        return np.zeros(1)
    
    # (1)  FFT transform
    y=np.copy(xin)
    # No windows .... ok
    windows0= signal.hamming(len(y))  # hanning windows
    yw= y # * windows0
    yf = fftpack.fft(yw)
    freqList = fftpack.fftfreq( len(y), d=1.0/ sr)
    # 
    yf_abs= np.abs(yf)   # get magnitude
    yf_phase= np.angle(yf) # get phase
    
    # (2) FFT transform
    ry=np.copy(ref_phase)
    # No windows .... ok.
    ryw= ry # * windows0
    ryf = fftpack.fft(ryw)
    # 
    ryf_abs= np.abs(ryf)   # get magnitude
    ryf_phase= np.angle(ryf) # get phase
    
    
    # (3) inverse FFT Transform
    yf2= yf_abs * np.exp(1.0j* ryf_phase)  # apply ref_phase's phase information
    y2 = np.real(fftpack.ifft(yf2))  # * windows0 if overlap-add with windows
    
    
    if show:
        # plot whole waveform
        fig = plt.figure()
        
        ax0 = fig.add_subplot(211)
        plt.xlabel('mSec')
        plt.ylabel('level')
        plt.title( 'waveform: blue original: red reference: yellow phase changed' )
        plt.plot( (np.arange(len(y)) * 1000.0 / sr) ,xin, color='b')
        plt.plot( (np.arange(len(ref_phase)) * 1000.0 / sr) ,ref_phase, color='r')
        plt.plot( (np.arange(len(y2)) * 1000.0 / sr) ,y2, color='y')
        plt.grid()
        
        ax1 = fig.add_subplot(212)
        plt.title('frequency response' + comment0)
        plt.ylabel('Amplitude [dB] solid line')
        plt.xlabel('Frequency [Hz]')
        plt.semilogx(freqList[1: int(len(yf)/2)], 20.0 * np.log10(np.abs(yf[1: int(len(yf)/2)])), color='b',)
        plt.semilogx(freqList[1: int(len(yf)/2)], 20.0 * np.log10(np.abs(ryf[1: int(len(yf)/2)])), color='r',)
        
        ax2 = ax1.twinx()
        angles = np.unwrap(yf_phase)
        angles = angles / ((2.0 * np.pi) / 360.0)
        rangles = np.unwrap(ryf_phase)
        rangles = rangles / ((2.0 * np.pi) / 360.0)
        plt.semilogx(freqList[1: int(len(yf)/2)], angles[1: int(len(yf)/2)], color='b', linestyle='dashed')  # plt.plot(flist, angles, 'g')
        plt.semilogx(freqList[1: int(len(yf)/2)], rangles[1: int(len(yf)/2)], color='r', linestyle='dashed')  # plt.plot(flist, angles, 'g')
        plt.ylabel('Angle(deg) dash line')
        
        
        plt.grid()
        fig.tight_layout()
        plt.show()
        
    return y2
