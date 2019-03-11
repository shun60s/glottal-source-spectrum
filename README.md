# glottal source spectrum   

A trial estimation of glottal source spectrum by anti-formant filter and inverse radiation filter.  

[github repository](https://github.com/shun60s/glottal-source-spectrum/)  

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

estimate of glottal source spectrum condition by anti-formant filter and inverse radiation filter  
```
examples:
  vowel /a/        : python3 est_gss1.py -w a_1-16k.wav  -f 3 -g 1  
  nasal voice /na/ : python3 est_gss1.py -w na_1-16k.wav -f 8 -g 2  
```
![figure4](docs/inverse_filter_output_a_3.png)  

glottal source spectrum condition  
![figure5](docs/source_frequecny_response_a_3.png)  


resampling wav to 16Khz sampling  
```
python resample1.py -w wav-file-name(mono,16bit)  
```


## License    
MIT  
except LPC.py  

