#coding:utf-8

# return candidate position set of one pitch duration near center of the frame
# by differential change point and threshold from bottom line.
# return 0 if there is no.
#
# 中心付近の1ピッチ分の候補インデックス[sp,ep]を返す。
# 候補が無いときは零を返す。
#
# 微分の変化点と閾値により候補を選出する。

import numpy as np
import matplotlib.pyplot as plt

# Check version
#  Python 3.6.4, 64bit on Win32 (Windows 10)
#  numpy (1.14.0)


def diff_ana(y, sr, show=False):
	# (1)　傾きの変化より選択
    f_prime=np.gradient(y) # 数値勾配（傾き）
    indices_diff0 = np.where( np.diff(np.sign(f_prime)) > 0.0 )[0]    # 符号(-1,0,1)化したものの差分をとり、正値の変化点を検出する
    # (2)　底辺に近い値を選択
    thres0= (np.amax(y) - np.amin(y)) * 0.25 + np.amin(y)  # 最小値から振幅幅の25％までの値を候補として使う。
    indices_thres0 = np.where( y < thres0 )[0]
    # (3)  上記の条件を満たす 論理積 を取る
    indices=np.sort(np.array(list( set(indices_diff0) & set(indices_thres0))))
    infections = y[indices]
    
    if len(indices) >= 2: #　候補が2個以上のときに、探す。
        index0= np.argmin(np.abs(indices - len(y)/2))  #　中心に一番近いインデックスを求める
        
        if len(indices) == 2: # 候補が2個しかないときは
            sp= indices[0]
            ep= indices[1]
        elif index0 < len(y)/2 and indices[-1] > len(y)/2 :  # そのインデックスが中心より前ならば
    	    sp= indices[index0]
    	    ep= indices[index0+1]
        else:
            sp= indices[index0-1]
            ep= indices[index0]
    else:  # 候補が無い
        sp=0
        ep=0
    
    indices1=np.array([sp,ep])
    infections1 = y[indices1]
    
    #print ( indices, indices1)
    #print ('select index, [Hz]', indices1, (sr / (indices1[1]-indices1[0])) )
    
    
    if show:
        fig = plt.figure()
        ax1 = fig.add_subplot(311)
        plt.title('diff: two red cirles shows selected portion')
        plt.xlabel('mSec')
        plt.ylabel('level')
        ax1.plot(np.arange(len(y)) * 1000.0 / sr, y, 'bo-', ms=2)
        ax1.plot(indices * 1000.0 / sr, infections, 'yo', ms=5)
        ax1.plot(indices1 * 1000.0 / sr, infections1, 'ro', ms=5)
        
        ax2 = fig.add_subplot(312)
        ax2.plot(np.arange(len(f_prime)) * 1000.0 / sr, f_prime, 'ro', ms=5)
        
        ax3 = fig.add_subplot(313)
        f_prime2=np.gradient(f_prime) 
        indices2 = np.where(np.diff(np.sign(f_prime2)))[0]
        infections2 = y[indices2]
        ax3.plot(np.arange(len(y)) * 1000.0 / sr, y, 'bo-', ms=2)
        ax3.plot(indices2 * 1000.0 / sr, infections2, 'ro', ms=5)
        
        plt.show()
    
    return int(sp), int(ep)
