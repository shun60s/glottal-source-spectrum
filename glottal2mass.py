#coding:utf-8

#-----------------------------------------------------------------------------------------------
# modified A.E.Rosenberg's formula to simply resemble two-mass model developed by Ishizaka and Flanagan
#
#  y=  max( yg( t + dt ), yg(t) * gain0)
#  dt is time difference.  If dt=0, then one pluse only
#  gain0 is adjustment factor of 2nd pulse amplitude to mix
#  Whole pulse length becomes trise + tfall + tdiff
#
#------------------------------------------------------------------------------------------------
# glottal voice source as input of Two Tubes Model of vocal tract
# Glottal Volume Velocity 
# based on A.E.Rosenberg's formula as Glottal Volume Velocity

import numpy as np
from scipy import signal
from matplotlib import pyplot as plt


# Check version
#  Python 3.6.4 on win32 (Windows 10)
#  numpy 1.14.0
#  scipy 1.0.0
#  matplotlib  2.1.1


class Class_Glottal2(object):
	def __init__(self, length_points= None, tclosed=5.0, trise=6.0, tfall=2.0, tdiff=0.0, gain0=1.0, sampling_rate=16000):
		# If length_points is specified, adjust real tclosed, trise, and tfall time, to meet whole length is length_points.
		# 比は入力のままで、長さがlength_pointsになるように実際のtclosed, trise, tfallの長さを調整する。
		#
		# initalize
		self.sr= sampling_rate
		if length_points is not None:
			self.length_points=length_points
			f=    ( length_points * 1000.0 / self.sr) / (tclosed + trise + tfall)
			self.tclosed=tclosed * f # duration time of close state [mSec]
			self.trise=trise * f     # duration time of opening [mSec]
			self.tfall=tfall * f     # duration time of closing [mSec]	
			self.tdiff=tdiff * f     # time difference of two pulse [mSec]			
			self.N1=int( (self.tclosed / 1000.) * self.sr )
			self.N2=int( (self.trise / 1000.) * self.sr )
			self.N3=int( (self.tfall / 1000.) * self.sr )
			self.dt=int( (self.tdiff / 1000.) * self.sr )
			self.LL= self.N1+ self.N2 + self.N3
			
			if self.LL > self.length_points or self.LL < (self.length_points -3):
				print ('error: (1)')
				print ('length_points - LL', self.length_points - self.LL)
			elif (self.length_points - self.LL) >= 1:
				self.N1 +=1
				self.tclosed= self.N1 * 1000 / self.sr
				if (self.length_points - self.LL) >= 2:
					self.N2 +=1
					self.trise= self.N2 * 1000 / self.sr
					if (self.length_points - self.LL) >= 3:
						self.N3 +=1
						self.tfall= self.N3 * 1000 / self.sr
			self.LL= self.N1+ self.N2 + self.N3
			# check
			if self.length_points !=  self.LL:
				print ('error: (2) ')
				print ('length_points - LL', self.length_points - self.LL)
		else:
			self.tclosed=tclosed  # duration time of close state [mSec]
			self.trise=trise      # duration time of opening [mSec]
			self.tfall=tfall      # duration time of closing [mSec]	
			self.tdiff=tdiff      # time difference of two pulse [mSec]
			self.N1=int( (self.tclosed / 1000.) * self.sr )
			self.N2=int( (self.trise / 1000.) * self.sr )
			self.N3=int( (self.tfall / 1000.) * self.sr )
			self.dt=int( (self.tdiff / 1000.) * self.sr )
			self.LL= self.N1+ self.N2 + self.N3
			self.length_points= self.LL
		
		self.gain0=gain0
		self.yg=self.make_one_plus()
		
	def make_one_plus(self,):
		# output yg
		yg=np.zeros(self.LL)
		yg1=np.zeros(self.LL)
		
		# at first one pulse generation
		for t0 in range(self.LL):
			if t0 < self.N1 :
				pass
			elif t0 <= (self.N2 + self.N1):
				yg1[t0]=  0.5 * ( 1.0 - np.cos( ( np.pi / self.N2 ) * (t0 - self.N1)) )
			else:
				yg1[t0]=  np.cos( ( np.pi / ( 2.0 * self.N3 )) * ( t0 - (self.N2 + self.N1) )  )
		
		# next product time difference pulse
		for t0 in range(self.LL):
			if int(t0 + self.dt) >= self.LL:
				yg[t0]=  yg1[t0] * self.gain0
			else:
				#yg[t0]= np.sqrt(yg1[t0]) * np.sqrt( yg1[int(t0 - self.dt)] )
				yg[t0]= np.amax( [ yg1[int(t0 + self.dt)], yg1[t0] * self.gain0 ])		
		
		# b=Numerator(bunsi), a=denominator(bunbo)
		self.a=np.ones(1)
		self.b=np.copy( yg[self.N1- self.dt :])
		self.b /= np.sum(self.b)  # in order to dc magnification is 1
		return yg

	def make_N_repeat(self, repeat_num=3):
		yg_repeat=np.zeros( len(self.yg) * repeat_num)
		for loop in range( repeat_num):
			yg_repeat[len(self.yg)*loop:len(self.yg)*(loop+1)]= self.yg
		return  yg_repeat
	
	def fone(self, f):
		# calculate one point of frequecny response
		xw= 2.0 * np.pi * f / self.sr
		yi=0.0
		yb=0.0
		for v in range (0,(self.dt + self.N2 + self.N3)):
			yi+=  self.yg[self.N1 - self.dt + v] * np.exp(-1j * xw * v)
			yb+=  self.yg[self.N1 - self.dt + v]
		val= yi/yb   # in order to dc magnification is 1
		return np.sqrt(val.real ** 2 + val.imag ** 2)
	
	def H0(self, freq_low=100, freq_high=5000, Band_num=256):
		# get Log scale frequecny response, from freq_low to freq_high, Band_num points
		self.amp=[]
		freq=[]
		self.bands= np.zeros(Band_num+1)
		fcl=freq_low * 1.0    # convert to float
		fch=freq_high * 1.0   # convert to float
		delta1=np.power(fch/fcl, 1.0 / (Band_num)) # Log Scale
		self.bands[0]=fcl
		#print ("i,band = 0", bands[0])
		for i in range(1, Band_num+1):
			self.bands[i]= self.bands[i-1] * delta1
			#print ("i,band =", i, bands[i]) 
		for f in self.bands:
			self.amp.append(self.fone(f) )
		
		return   np.log10(self.amp) * 20, self.bands # = amp value, freq list
		
	def f_show(self, worN=1024):
		# show frequency response, using scipy
		wlist, fres = signal.freqz(self.b, self.a, worN=worN)
		
		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		flist = wlist / ((2.0 * np.pi) / self.sr)
		plt.title('frequency response')
		ax1 = fig.add_subplot(111)
		
		plt.semilogx(flist, 20 * np.log10(abs(fres)), 'b')  # plt.plot(flist, 20 * np.log10(abs(fres)), 'b')
		
		plt.semilogx(self.bands, 20 * np.log10(np.array(self.amp)), 'r')  # plt.plot(flist, 20 * np.log10(abs(fres)), 'b')
		
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


if __name__ == '__main__':
	
	#  parameters variation list 
	tclosed_list = [ 5.0, 5.0 , 6.0, 5.0, 7.0, 7.0]
	trise_list   = [ 6.0, 6.0 , 5.0, 6.0, 5.0, 5.0]
	tfall_list   = [ 2.0, 2.0 , 2.0, 2.0, 1.0, 1.0]
	tdiff_list   = [ 0.0, 1.0 , 1.0, 2.0, 0.0, 1.0]
	gain0_list   = [ 1.0, 1.1 , 1.1, 1.1, 1.0, 1.05]
	
	# draw waveform and its frequency response
	color_list = ["r", "g", "b", "c", "m", "y", "k", "w"]
	fig = plt.figure()
	
	for i in range( len(tclosed_list) ):
		
		# instance
		glo=Class_Glottal2(tclosed=tclosed_list[i], trise=trise_list[i], tfall=tfall_list[i], tdiff=tdiff_list[i], gain0=gain0_list[i])
		wlist, fres = signal.freqz(glo.b, glo.a, worN=1024)
		flist = wlist / ((2.0 * np.pi) / glo.sr)
		
		ax1 = fig.add_subplot(211)
		plt.xlabel('mSec')
		plt.ylabel('level')
		plt.title("Glottal waveform by modified A.E.Rosenberg's formula")
		plt.plot( (np.arange(len(glo.yg)) * 1000.0 / glo.sr) , glo.yg, color=color_list[i])
		#plt.grid()
		
		ax1 = fig.add_subplot(212)
		plt.title('frequency response')
		plt.ylabel('Amplitude [dB]')
		plt.xlabel('Frequency [Hz]')
		plt.semilogx(flist, 20 * np.log10(abs(fres)), color=color_list[i])  # plt.plot(flist, 20 * np.log10(abs(fres)),color=color_list[i])
		#plt.grid()
		
	ax1 = fig.add_subplot(211)
	plt.grid()	
	ax1 = fig.add_subplot(212)
	plt.grid()	
	fig.tight_layout()
	plt.show()
	
	
	
#This file uses TAB
