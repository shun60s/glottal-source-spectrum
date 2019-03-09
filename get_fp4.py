#coding: utf-8

# estimate formant peak frequency , Q factor, and pitch frequency based on LPC analysis
# return candidates of format peak frequecny, Q factor, and Log-scale LPC spectral
# input wav must be mono, 16bit.
#
# LPC分析のよるホルマント周波数とピッチ周波数の推定。
# 対数LPCスペクトとホルマント周波数の候補とＱの推定値を返す。
# 入力のWAVファイルは16ビットのMONO信号を想定している。
#
# 2018/4/20 出力形式を一部変更　suffix 2へ変更
# 2019/2/22 出力形式を一部変更　suffix 3へ変更
# 2019/3/6  出力形式を一部変更　suffix 4へ変更

import wave
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt

from LPC import *

# Check version
#  Python 3.6.4, 64bit on Win32 (Windows 10)
#  numpy (1.14.0)
#  matplotlib  2.1.1
#  scipy (1.0.0)

class Class_get_fp(object):
	def __init__(self,NFRAME=640, NSHIFT=320, lpcOrder=32, Delta_Freq=5, max_num_formants=5):
		
		self.NFRAME=NFRAME     # 640 sr=16Khz 40mS  # 400 sr=16Khz 25mS 
		self.NSHIFT=NSHIFT     # 320 sr=16Khz 20mS  # 160 sr=16khz 10mS
		self.lpcOrder=lpcOrder
		self.Delta_Freq = Delta_Freq  # frequency resolution to compute frequency response
		self.window = np.hamming(self.NFRAME)  # Windows is Hamming
		self.preemph=0.97  # pre-emphasis 
		self.max_num_formants =max_num_formants   # maximum number of formant candidate to detect


	def get_fp(self,file_name, frame_num=None , figure_show=False ):
		#   入力：wave ファイル　mono 16bit
		#
		#   出力：LPC対数スペクト周波数の行列
		#         ホルマント周波数の候補のインデックス
		#         ピッチ周波数の候補のインデックス
		#         
		#   もし、frame_numが指定されていた場合は、そのフレーム番号だけ計算する
		
    	# read wave file
		waveFile= wave.open( file_name, 'r')
		
		nchannles= waveFile.getnchannels()
		samplewidth = waveFile.getsampwidth()
		sampling_rate = waveFile.getframerate()
		nframes = waveFile.getnframes()
		self.FreqPoints= int((sampling_rate / 2.0) / self.Delta_Freq)
		self.df0 = (sampling_rate /2.) / self.FreqPoints
		self.dt0 = 1.0 / sampling_rate
		
		# check input wave file condition
		assert nchannles == 1, ' channel is not MONO '
		assert samplewidth==2, ' sample width is not 16bit '
		
		
		if frame_num is not None and frame_num >= 0:
			waveFile.readframes( self.NSHIFT * frame_num ) # dummy read
			nframes=self.NFRAME
			buf = waveFile.readframes(nframes ) # read only frame_num portion
			self.frame_num=frame_num
		else:
			buf = waveFile.readframes(-1) # read all at oance
			self.frame_num=None
		
		waveFile.close()
		
		# 16bit integer to float32
		data = np.frombuffer(buf, dtype='int16')
		fdata = data.astype(np.float32)
		
		count= int(((nframes - ( self.NFRAME - self.NSHIFT)) / self.NSHIFT))
		
		
		# prepare output
		spec_out= np.zeros([count,self.FreqPoints])
		fout = np.zeros([count,self.max_num_formants])
		fout_index = np.ones([count,self.max_num_formants]) * -1
		pout = np.zeros(count)
		pout_index = np.ones(count) * -1
		
		Qout= np.zeros([count,self.max_num_formants])
		Low_index= np.zeros([count,self.max_num_formants])
		High_index=np.zeros([count,self.max_num_formants])
		
		pos = 0  # position
		countr=0
		
		for loop in range(count):
			
			## copy to avoid original over-change
			frame = fdata[pos:pos + self.NFRAME].copy()
			
			## pre-emphasis
			frame -= np.hstack((frame[0], frame[:-1])) * self.preemph
			## do window
			windowed = self.window * frame
			## get lpc coefficients
			a,e=lpc(windowed, self.lpcOrder)
			## get lpc spectrum
			w, h = scipy.signal.freqz(np.sqrt(e), a, self.FreqPoints)  # from 0 to the Nyquist frequency
			lpcspec = np.abs(h)
			lpcspec[lpcspec < 1.0] = 1.0  # to avoid log(0) error
			loglpcspec = 20 * np.log10(lpcspec)
			spec_out[loop]=loglpcspec # store to output
			## get formant candidate
			f_result, i_result=self.formant_detect(loglpcspec, self.df0)
			## sort into max_num_formants
			if len(f_result) > self.max_num_formants:
				fout[loop]=f_result[0:self.max_num_formants]
				fout_index[loop]=i_result[0:self.max_num_formants]
			else:
				fout[loop]=f_result[0:len(f_result)]
				fout_index[loop]=i_result[0:len(f_result)]
			
			## get Q factor
			Q_list,_,_, Low_index0, High_index0= self.bandwidth_detect( loglpcspec, self.df0, fout_index[loop])
			Qout[loop]=Q_list
			Low_index[loop]= Low_index0
			High_index[loop]= High_index0
			
			if figure_show:
				self.show_Q_point(spec_out[loop], fout_index[loop], Low_index[loop], High_index[loop], self.df0, loop)
			
			## calcuate lpc residual error (= input source)
			r_err=residual_error(a, windowed)
			## autocorrelation of lpc residual error (= input source)
			a_r_err=autocorr(r_err)
			a_f_result, a_i_result = self.pitch_detect(a_r_err, self.dt0)
			if len(a_f_result) > 0: # if candidate exist,
				pout[loop]=a_f_result[0]
				pout_index[loop]=a_i_result[0]
			
			
			## print output of candidates of [formants],  frequency[Hz]
			"""
			if countr == 0:
				print ('candidates of [formants],  frequency[Hz] ')
			print (fout[loop])
			"""
			
			# index count up
			countr +=1
			# next
			pos += self.NSHIFT
		
		return spec_out, fout, pout, Qout, fout_index, Low_index, High_index


	def formant_detect(self,input0, df0, f_min=250, db_min=3):
		#   対数スペクトルから
		#   山型（凸）のピークポイントを見つける
		#
		#   入力：対数スペクトル
		#         周波数単位
		#         （オプション）最低の周波数
		#         （オプション）周囲からdb_min dB以上高い場合をピークとみなす
		#
		#   出力：ピークのインデックス
		#         ピークの周波数
		is_find_first= False
		f_result=[]
		i_result=[]
		for i in range (1,len(input0)-1):
			if f_min is not None and  df0 * i <= f_min :
				continue
			if input0[i] > input0[i-1] and input0[i] > input0[i+1] :
				if not is_find_first :
					f_result.append( df0 * i)
					i_result.append(i)
					is_find_first =True
				else:
					f_result.append( df0 * i)
					i_result.append(i)
		
		# 周囲からdb_min dB以上高い場合をピークとみなす
		Q_list,_,_,_,_ =  self.bandwidth_detect(input0, df0, i_result, db_min=3)
		f_result2= np.array(f_result)[np.where( np.array(Q_list) > 0.0)[0]]
		i_result2= np.array(i_result)[np.where( np.array(Q_list) > 0.0)[0]]
		
		return list(f_result2), list(i_result2)  # 旧品との互換性を保つためリストに戻す
		
		
	def bandwidth_detect(self,input0, df0, peak_index_list, db_min=3):
		#   対数スペクトルのピークw0から
		#   3dB(db_minのデフォルト値)低下した周波数を求める
		#
		#   入力：対数スペクトル  * 高精度の求めるには　周波数分解能が高いことが求められる。
		#         周波数単位
		#         ピークのインデックス
		#
		#   出力：-3dB低下した周波数w1（ピークより低い周波数の方）
		#         -3dB低下した周波数w2（ピークより高い周波数の方）
		#          Q= w0/(w2-w1) 定義できないときは零が入る
		#          暫定的に　w1 又は w2のどちらかが分かっているときは、Q= w0/(2*(w0-wx))を入れておく
		Q_list=[]
		Low_freq_list=[]
		High_freq_list=[]
		Low_freq_index=[]
		High_freq_index=[]
		
		for ipk in (peak_index_list):
			low_index= 0
			high_index= 0
			f_peak= input0[int(ipk)]
			f_peak_3dB= input0[int(ipk)] - db_min  # 1/sqrt(2) => -3dB
			
			# ピークより低い周波数の方を探す
			for i in range (int(ipk),0,-1):
				if input0[i] <= f_peak_3dB:
					low_index= i
					break
				elif input0[i] > f_peak:
					break
					
			# ピークより高い周波数の方を探す
			for i in range (int(ipk),len(input0)):
				if input0[i] <= f_peak_3dB:
					high_index= i
					break
				elif input0[i] > f_peak:
					break
									
			#Q= w0/(w2-w1)の計算　 定義できないときは零が入る
			if low_index > 0 and high_index > 0:
				Q= 1.0 * int(ipk) / (high_index - low_index)
			elif low_index > 0:
				Q= 1.0 * int(ipk) / (2.0 * (int(ipk) - low_index))  # 暫定値を入れておく
			elif high_index > 0:
				Q= 1.0 * int(ipk) / (2.0 * ( high_index - int(ipk)))  # 暫定値を入れておく
			else:
				Q=0.0
		
			Q_list.append( Q )
			Low_freq_list.append( df0 * low_index)
			High_freq_list.append( df0 * high_index)
			Low_freq_index.append( low_index)
			High_freq_index.append( high_index)
			
			#print ( Q, df0 * low_index, df0 * high_index, df0 * int(ipk))
		
		return Q_list, Low_freq_list, High_freq_list, Low_freq_index, High_freq_index
		
		
	def pitch_detect(self, input0, dt0, ratio0=0.3, f_min=100, f_max=500):
		# 　自己相関の
		# 　山と谷の両方のピークを求める
		#
		#   入力：lpc予測残差の自己相関
		#         時間単位
		#         （オプション）自己エネルギー0次成分に対する比率（これ以上を対象とする）
		#         （オプション）最低の周波数
		#         （オプション）最大の周波数
		#
		#   出力：最大ピークのインデックス
		#         最大ピークの周波数の値
		#
		#
		is_find_first= False
		f_result=[]
		i_result=[]
		v_result=[]
		for i in range (1,len(input0)-1):
			if np.abs(input0[i]) < np.abs(input0[0] * ratio0):
				continue
			fp= 1.0 / (dt0 * i)
			if f_max is not None  and fp >= f_max :
				continue
			if f_min is not None and  fp <= f_min :
				continue
			if input0[i] > input0[i-1] and input0[i] > input0[i+1] :
				if not is_find_first :
					f_result.append( fp)
					i_result.append(i)
					v_result.append( input0[i])
					is_find_first =True
				else:
					f_result.append( fp)
					i_result.append(i)
					v_result.append( input0[i])
			elif input0[i] < input0[i-1] and input0[i] < input0[i+1] :
				if not is_find_first :
					f_result.append( fp)
					i_result.append(i)
					v_result.append( input0[i] )
					is_find_first =True
				else:
					f_result.append( fp)
					i_result.append(i)
					v_result.append( input0[i])
		
		if is_find_first:  # 最大のピークを探す
			a=np.argmax( np.array(v_result))
			f_result2= [ f_result[np.argmax( np.array(v_result))] ]
			i_result2= [ i_result[np.argmax( np.array(v_result))] ]
		else: #　候補なし
			f_result2=[]
			i_result2=[]
			
		return f_result2, i_result2
	
	def show_Q_point(self,spec_out, fout_index_list, low_index_list, high_index_list, df0, loop):
		# draw peak, -3dB points on frequency response
		fout_index=np.array(fout_index_list, dtype=np.int)
		low_index=np.array(low_index_list, dtype=np.int)
		high_index=np.array(high_index_list, dtype=np.int)
		#
		fig = plt.figure()
		ax1 = fig.add_subplot(211)
		if self.frame_num is None:
			index0=loop
		else:
			index0=self.frame_num
		plt.title('frame no. ' + str(index0) + ': frequency response')
		plt.xlabel('Frequency [Hz]')
		plt.ylabel('Amplitude [dB]')
		ax1.plot(np.arange(len(spec_out)) * df0, spec_out, 'b', ms=2)
		ax1.plot( fout_index * df0 , spec_out[fout_index], 'ro', ms=3)
		ax1.plot( low_index[ low_index > 0] * df0, spec_out[low_index[ low_index > 0]], 'yo', ms=3)
		ax1.plot( high_index[ high_index > 0] * df0, spec_out[high_index[ high_index > 0]], 'yo', ms=3)
		plt.grid()
		plt.axis('tight')
		plt.show()


if __name__ == '__main__':
	
	import argparse
	
	parser = argparse.ArgumentParser(description='estimation formant peak frequency and Q factor')
	parser.add_argument('--wav_file', '-w', default='a_1-16k.wav', help='specify a wav-file-name(mono,16bit)')
	parser.add_argument('--frame', '-f', type=int, default=-1, help='specify the frame number, set negative value if ignore')
	args = parser.parse_args()
	
	# instance
	fp0=Class_get_fp()
	
	# get
	spec_out, fout, pout, Qout, fout_index, Low_index, High_index = fp0.get_fp( args.wav_file, frame_num=args.frame , figure_show=True )
	
	print ( fout)  # candidates of format peak frequecny
	print ( Qout)  # Q factor of them
	
	
#This file uses TAB
