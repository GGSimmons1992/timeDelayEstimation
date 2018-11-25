#Used for testing new codes/algorithims
import numpy as np
import random
import matplotlib.pyplot as plt
import TDEBasic as TDE
"""
#Test 1: x in range(0,len(array))
c=[0,1,2,3,4,5,6,7,8,9,10]
for x in range(0,len(c)):
    print x
print "\nlength of c is {}".format(len(c))
"""

"""
#Test 2a: Splicing and matching
a=[]
for x in range(0,11):
    a.append(x)
a=np.array(a)
starter=5
block=3
shift=2
pad=starter%shift
print pad
CorrSize=int(starter/shift)+int((len(a)-(starter+block))/shift)
print CorrSize
miniStart=pad
for x in range(0,CorrSize+1):
    b=a[miniStart:miniStart+block]
    print "b{}\n{}\n".format(x,b)
    miniStart=miniStart+shift
#Test 2b: finding dt
"""
"""
Assume correlation is found in b1
Then dt should be -0.5+/-1.5
"""
"""
i=1
errorDT=block/2.0
dt=pad+(i*shift)+errorDT-starter
print "dt is: {}+/-{}".format(dt,errorDT)
"""
"""
#Test 3) shifting arrays
def Correlation(miniA,miniB):
    miniA=np.array(miniA)-np.mean(miniA)
    miniB=np.array(miniB)-np.mean(miniB)
    Corre=np.sum(miniA*miniB)/(len(miniA)*np.std(miniA)*np.std(miniB))
    #print "{}, {}, {}".format(miniA,miniB,Corre)
    return Corre


t=[]
A=[]
B=[]
C=[]
bAccuracy=[0.0]*10
cAccuracy=[0.0]*10
for y in range(0,2000):
        t.append(y)
        A.append(0)
        B.append(0)
        C.append(0)
t=np.array(t)
A=np.array(A)
B=np.array(B)
C=np.array(C)
start=1000
block=400
for x in range(0,10):
    Amp=(random.random()*100)-(random.random()*100)
    f=random.random()*99+1
    shift=int((random.random()*500)+0.4999)-int((random.random()*500)+0.4999)
    for z in range(0,len(t)):
        A[z]=(Amp*np.sin(f*z))
        B[z]=(Amp*np.sin(f*z+shift))
        C[z]=(Amp*np.sin(f*(z+shift)))
"""
""" 
    A=np.array(A)-np.mean(A)
    B=np.array(B)-np.mean(B)
    C=np.array(C)-np.mean(C)
"""
"""
    miniA=A[start:start+block]
    miniB=B[start+shift:start+block+shift]
    miniC=C[start+shift:start+block+shift]
    #print "\n{}sin({}t) and {}sin({}t+{})".format(Amp,f,Amp,f,shift)
    bAccuracy[x]=Correlation(miniA,miniB)
    #print "{}sin({}t) and {}sin({}(t+{}))".format(Amp,f,Amp,f,shift)
    cAccuracy[x]=Correlation(miniA,miniC)
print "\n"
print bAccuracy
print cAccuracy
print "\nMethod B (f*t+shift) is {}% accurate".format(np.mean(bAccuracy)*100)
print "\nMethod C (f*(t+shift)) is {}% accurate".format(np.mean(cAccuracy)*100) 
"""
#Test 4) plotting wavefunctions
sineA,sineB,t,tShift,Amp,omega,shift=TDE.shiftedSines()
length=len(sineA)
noiseA,noiseB,SNR=TDE.NoiseMaker(length,Amp)
print 'SNR is {}'.format(SNR)
packetA,packetB,sigma=TDE.Packets(t,tShift,length)
wavepacketA=(sineA*packetA)+noiseA
wavepacketB=(sineB*packetB)+noiseB
plt.subplot(3,1,1)
plt.plot(t,sineA,'b',t,noiseA,'r',t,packetA,'g')
plt.subplot(3,1,2)
plt.plot(t,sineA,'b',t,noiseA,'r',t,packetA,'g')
plt.subplot(3,1,3)
plt.plot(t,wavepacketA,'b',t,wavepacketB,'r')
plt.show()
