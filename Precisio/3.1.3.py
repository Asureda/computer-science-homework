# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 10:59:36 2019

@author: orica
"""
#import wave,struct
from scipy.io import wavfile
#import soundfile as sf
#w_in=wave.open("ona.wav","r")
#w_in=open("ona.wav","r")
f_out=open("ona_out.dat","w")
f_in=open("ona_in.dat","w")
t=[]
w=[]
#leng=w_in.getnframes()
#for line in range(0,leng-2):
sr,data=wavfile.read('ona.wav')
#data,sr=sf.read('ona.wav')
print(data,sr)
for line in w_in:
    #l=w_in.readframes(1)
    #data=struct.unpack("<h",l)
    fs,data=w_in.read()
    print(data)
    t.apend(line[0])
    w.apend(line[1])
for i in range(len(w)-2):
    f_in.write(str(t[i])+'\t'+str(w[i]))
    if ((i%2) != 0):
        aux=(w[i-1]+w[i+1])/2.0
        f_out.write(str(t[i])+'\t'+str(aux))
