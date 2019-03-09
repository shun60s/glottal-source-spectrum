# glottal source spectrum   



## content  

some frequency response  samples of pseudo glottal source waveform  
```
python3 glottal2mass.py
```
![figure1](docs/glottalwaves2freqres.png) 


get a low shelving filter as inverse filter against high pass filter that simulates radiation from mouth  
```
python3 iir_LowShelving1.py
```
![figure2](docs/lowShelving2HPF.png)  



estimate formant peak frequency and Q factor based on LPC analysis  
```
python3 get_fp4.py -w wav-file-name(mono,16bit) -f frame-number
```
![figure3](docs/formant_and_Q-3dB_points.png)  


resampling to 16Khz sampling  
```
python resample1.py -w wav-file-name(mono,16bit)  
```


## License    
MIT  
except LPC.py  

