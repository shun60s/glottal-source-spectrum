#coding:utf-8

# trial estimation of glottal source spectrum condition by anti-formant filter and inverse radiation filter 
# under following hypotheses.
#  (1) glottal source spectrum characterizes descending rightwards without sharp peak.
#  (2) resonance strength of formant is roughly same regardless of formant frequency.

#　声門の音源スペクトルは　険しいピークがない右肩下がりの特性である
#  フォルマントの共鳴の強さは周波数によらず同じ程度である
#  の仮説のもと、
#  口の放射特性の逆フィルターとフォルマント周波数で減衰するフィルターを使って
#  声門の音源のスペクトルの状態を予想する

import argparse
import numpy as np
from scipy import signal
from scipy.io.wavfile import read as wavread
from matplotlib import pyplot as plt

from get_fp4 import *
from BPF import *
from iir_LowShelving1 import *
from iir_peaking1     import *
from diff_ana         import *
from fft_ana          import *
from glottal2mass     import *

# Check version
#  Python 3.6.4 on win32 (Windows 10)
#  numpy 1.14.0 
#  matplotlib  2.1.1
#  scipy 1.0.0

class Class_estimate_gss1(object):
    def __init__(self, path0):  # sampling_rate=16000):
        # initalize
        sr, y = wavread(path0)
        self.yg= y / (2 ** 15)
        self.sr= sr
        print ('sampling rate ', sr)
        
        self.fp0=Class_get_fp()
        # calculate lpc log-spectrum and  formant, Q 
        spec_out, fout, pout, Qout, fout_index, Low_index, High_index = self.fp0.get_fp(path0)
        self.fout= fout
        self.Qout= Qout
        self.gain= np.ones( self.Qout.shape)
        self.NFRAME=self.fp0.NFRAME
        self.NSHIFT=self.fp0.NSHIFT
        print ('formant freq and Q, per frame')
        print ( self.fout )
        print ( self.Qout)
        
        # apply a low shelving filter as inverse filter against high pass filter that simulates radiation from mouth
        self.iir_LS1=Class_IIR_LowShelving1(sampling_rate=self.sr)
        self.invhpf_wav= self.iir_LS1.filtering(self.yg)
        self.invhpf_wav /=  (np.amax(np.abs(self.invhpf_wav)) / np.amax(np.abs(self.yg))) # normalize
        
    def analysis(self, frame_num=None, gain_pattern=1, figure_show = True):
        
        # only process one frame, if frame_num is specified
        for l in range(self.fout.shape[0]):
            if frame_num is not None:
                if frame_num >= 0 and l != frame_num:
                    continue
            
            # compute start point and end point of current l-th frame
            sp= self.NSHIFT * l
            ep= sp + self.NFRAME
            if ep > len(self.yg):
                ep= len(self.yg)
            
            print ('frame no.', l, '  start[ms]', int(sp * 1000 / self.sr))
            # process BPF
            self.bpf1=Class_BPF(fc=self.fout[l,0], Q=self.Qout[l,0], sampling_rate=self.sr)
            self.bpf2=Class_BPF(fc=self.fout[l,1], Q=self.Qout[l,1], sampling_rate=self.sr)
            self.bpf3=Class_BPF(fc=self.fout[l,2], Q=self.Qout[l,2], sampling_rate=self.sr)
            self.bpf4=Class_BPF(fc=self.fout[l,3], Q=self.Qout[l,3], sampling_rate=self.sr)
            
            self.bpf_list=[self.bpf1, self.bpf2, self.bpf3, self.bpf4]
            
            # process BPF, filtering independently
            self.f1_wav=self.bpf_list[0].iir2(self.yg)
            self.f2_wav=self.bpf_list[1].iir2(self.yg)
            self.f3_wav=self.bpf_list[2].iir2(self.yg)
            self.f4_wav=self.bpf_list[3].iir2(self.yg)
            self.filtering_list=[ self.f1_wav, self.f2_wav, self.f3_wav, self.f4_wav]
            
            # set drop gain of iir peaking filter as anti-formant boost.
            # Try several assumed gain patterns to study appropriate drop gain.
            self.gain_pattern= gain_pattern
            if self.gain_pattern == 1:
                #　（パターン１）フォルマントの強さ(ゲイン)は、周波数によらず同じ程度と仮定する。
                #   (Pattern 1)resonance strength of formant is roughly the same regardless of the frequency.
                #  
                # all -20dB(=0.1) drop
                self.gain= np.ones( self.Qout.shape) * 0.1
                self.analysis_sub(l, sp, ep,  figure_show=figure_show)
            if self.gain_pattern == 2:
                #　（パターン２）鼻音効果で 2kHz以上のフォルマントの強さ（ゲイン）は弱まっている仮定とする。
                #  　　　　 2KHzは可変値。
                #   (Pattern 2)due to nose effect, resonance strength of the formant over 2kHz (adjustable) become weak.
                #   When formant frequency > 2KHz -10dB(=0.3162), other -20dB(=0.1) drop
                self.gain= np.ones( self.Qout.shape) * 0.1
                highside_gain=-10
                self.gain[np.where(self.fout > 2000)]=  np.power(10.0, ( highside_gain /20))
                self.analysis_sub(l, sp, ep, figure_show=figure_show)
            
    def analysis_sub(self, l, sp, ep, figure_show):
        # instance peaking drop filter
        self.pk1=Class_IIR_Peaking1(fpeak=self.fout[l,0], gain=self.gain[l,0], Q=self.Qout[l,0] , sampling_rate=self.sr)
        self.pk2=Class_IIR_Peaking1(fpeak=self.fout[l,1], gain=self.gain[l,1], Q=self.Qout[l,1] , sampling_rate=self.sr)
        self.pk3=Class_IIR_Peaking1(fpeak=self.fout[l,2], gain=self.gain[l,2], Q=self.Qout[l,2] , sampling_rate=self.sr)
        self.pk4=Class_IIR_Peaking1(fpeak=self.fout[l,3], gain=self.gain[l,3], Q=self.Qout[l,3] , sampling_rate=self.sr)
        self.pk_list=[ self.pk1, self.pk2, self.pk3, self.pk4]
        
        # process filtering in series
        self.pk1_wav=self.pk_list[0].filtering(self.invhpf_wav)
        self.pk2_wav=self.pk_list[1].filtering(self.pk1_wav)
        self.pk3_wav=self.pk_list[2].filtering(self.pk2_wav)
        self.pk4_wav=self.pk_list[3].filtering(self.pk3_wav)
        self.pk_filtering_list=[ self.pk1_wav, self.pk2_wav, self.pk3_wav, self.pk4_wav]
        
        # get a pitch duration in the frame, in order to avoid fundamental frequency F0 influence of frequency response
        # 基本周波数 F0の影響を除くため1ピッチ分の信号を取り出す
        sub_sp, sub_ep = diff_ana(   self.pk_filtering_list[-1] [sp:ep] , self.sr)
        
        #  show waveform
        if figure_show:
            self.plot_waveform2( l, sp, ep, sub_sp, sub_ep)
        
        # generate same length pseudo　glottal waveform as a reference
        # 同じ長さのリファレンス波形（基準とする波形）として、疑似的な声門の波形を生成する
        glo0=Class_Glottal2(length_points=(sub_ep-sub_sp), sampling_rate=self.sr)
        #glo0=Class_Glottal2(length_points=(sub_ep-sub_sp),tclosed=5.0, trise=5.0, tfall=0.8, tdiff=1.0, gain0=0.9, sampling_rate=self.sr)
        
        # get frequency response as glottal source spectrum
        comment0=': gain pattern  ' + str(self.gain_pattern)
        fft_ana( self.pk_filtering_list[-1][sp+sub_sp: sp+sub_ep], glo0.yg, self.sr, comment0, show=True)
        
        
    def plot_waveform2(self,loop, sp, ep, sub_sp, sub_ep):
        # plot every waveform per frame
        # set draw number
        max_display=2  # no display of BPF output
        fig = plt.figure()         
        
        """
        max_display=6  # display BPF output
        fig = plt.figure(figsize=(6, 7))  # adjust draw display area size
        """
        plt.subplot(max_display,1,1)
        plt.xlabel('mSec')
        plt.ylabel('level')
        plt.title( 'frame no. ' + str(loop) + ': blue original: red inverse radiation filter' )
        plt.plot( (np.arange(len(self.yg[sp:ep])) * 1000.0 / self.sr) , self.yg[sp:ep])
        plt.plot( (np.arange(len(self.invhpf_wav[sp:ep])) * 1000.0 / self.sr) , self.invhpf_wav[sp:ep], color='r')
        
        plt.subplot(max_display,1,2)
        plt.xlabel('mSec')
        plt.ylabel('level')
        plt.title( 'frame no. ' + str(loop) + ': anti-formant filter: red cirles, selected pitch portion' )
        plt.plot( (np.arange(len(self.pk_filtering_list[-1][sp:ep])) * 1000.0 / self.sr) , self.pk_filtering_list[-1][sp:ep], color='g')
        
        if sub_sp != 0 and sub_ep != 0 :
            indices1=np.array([sub_sp, sub_ep])
            infections1 = self.pk_filtering_list[-1][sp:ep][indices1]
            plt.plot( ( indices1 * 1000.0 / self.sr) , infections1, 'ro', ms=5)
        
        for i in range (max_display-2):
            plt.subplot(max_display,1,i+3)
            plt.xlabel('mSec')
            plt.ylabel('level')
            plt.title( 'f' + str(i+1) + ': '+ str(self.fout[loop,i]) + '[Hz] bpf output')
            plt.plot( (np.arange(len( self.filtering_list[i][sp:ep])) * 1000.0 / self.sr) ,  self.filtering_list[i][sp:ep])
            plt.grid()
            
        fig.tight_layout()
        plt.show()


if __name__ == '__main__':
    #
    parser = argparse.ArgumentParser(description='estimation glottal source spectrum condition')
    parser.add_argument('--wav_file', '-w', default='a_1-16k.wav', help='wav-file-name(mono,16bit)')
    parser.add_argument('--frame', '-f', type=int, default=-1, help='specify the frame number, set negative value if ignore')
    parser.add_argument('--gain', '-g', type=int, default=1, help='specify anti-formant drop gain pattern, set 1 if equal ')
    args = parser.parse_args()
    
    # examples:
    #  vowel /a/        : python3 est_gss1.py -w a_1-16k.wav  -f 3 -g 1
    #  nasal voice /na/ : python3 est_gss1.py -w na_1-16k.wav -f 8 -g 2
    
    # instance
    ana=Class_estimate_gss1(args.wav_file )
    ana.analysis( frame_num=args.frame, gain_pattern=args.gain )
    
    