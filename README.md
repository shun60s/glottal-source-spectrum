# glottal source spectrum   

A trial estimation of glottal source spectrum by anti-formant filter and inverse radiation filter.  

## contents  

some frequency response  samples of pseudo glottal source waveform  
```
python3 glottal2mass.py
```
![figure1](docs/glottalwaves2freqres.png) 


estimate formant peak frequency and Q factor based on LPC analysis  
```
python3 get_fp4.py -w wav-file-name(mono,16bit) -f frame-number
```
![figure2](docs/formant_and_Q-3dB_points-a-3.png)  


get a low shelving filter as inverse filter against high pass filter that simulates radiation from mouth  
```
python3 iir_LowShelving1.py
```
![figure3](docs/lowShelving2HPF.png)  


peaking filter class to drop formant boost portion, anti-formant filter  
```
iir_peaking1.py
```


resampling wav to 16Khz sampling  
```
python resample1.py -w wav-file-name(mono,16bit)  
```


## License    
MIT  
except LPC.py  

