"""
TDEBasic 2 point measurement
Used as backbone for multi-TDE calculation
"""
import numpy as np
import random
import matplotlib.pyplot as plt
import pylab
Debugger=0 #Assign non-zero value to turn on debugging checks noted by if(Debugger!=0).
pi=np.pi

def TDECorrelate(bigA,bigB,TimeStart,TimeBlock,TimeDelay,plotSwitch,theoryDT):
    #1)Initialize
    if(Debugger!=0):
        print 'TDEBasic.TDECorrelate:'
    bigA=np.array(bigA)-np.mean(bigA)
    bigB=np.array(bigB)-np.mean(bigB)
    NegCorrCells=int(TimeStart/TimeDelay)
    PosCorrCells=int((len(bigB)-(TimeStart+TimeBlock))/TimeDelay)
    TotalCorrCells=NegCorrCells+PosCorrCells
    Correlations=np.array([0.0]*TotalCorrCells)
    TimeShiftArray=np.array([0.0]*len(Correlations))
    pad=(TimeStart%TimeDelay)
    if(Debugger!=0):
        print "\npad is: {}\n".format(pad)
    miniA=bigA[TimeStart:(TimeStart+TimeBlock)]
    miniA=miniA-np.mean(miniA)
    starter=pad
    if(Debugger!=0):
        print len(Correlations)
    
    #2)DO ALL CORRELATIONS
    for x in range(0,len(Correlations)):
        TimeShiftArray[x]=starter-TimeStart
        if TimeShiftArray[x]==0:
            zero=x
        miniB=bigB[starter:(starter+TimeBlock)]
        miniB=miniB-np.mean(miniB)
        modulus=np.sum(miniA*miniB)
        CorreDenominator=(TimeBlock*np.std(miniA)*np.std(miniB))
        Correlations[x]=modulus/CorreDenominator
        del miniB
        starter=starter+TimeDelay
    
    #3)Analysis and Calculations
    TimeShiftArray=np.array(TimeShiftArray)
    i=(np.argmax(Correlations))
    if (Debugger!=0):
        print "\ni is:{}.\nmax coorelation is:{}\n".format(i,max(Correlations))
    dt=TimeShiftArray[i]
    percentCorrelation=max(Correlations)*100
    
    #Type 0 Error
    #errorDT=np.sqrt(TimeBlock**2+TimeDelay**2)
    
    #Type 1 Error
    errorDT=np.sqrt(abs(sum(np.square(Correlations*TimeShiftArray))/len(Correlations)))
    
    #Type 2 Error
    #absCorrelations=abs(Correlations)
    #errorDT=np.sqrt((sum(absCorrelations*np.square(TimeShiftArray)))/len(TimeShiftArray))
    
    #Type 3 Error
    """
    expectedT=sum(np.square(Correlations)*TimeShiftArray)/sum(np.square(Correlations))
    expectedTSquared=sum(np.square(Correlations*TimeShiftArray))/sum(np.square(Correlations))
    errorDT=np.sqrt(abs(expectedTSquared-np.square(expectedT)))
    """

    #Type 1a Error (Similar to 1, but it uses sum(Correlations) in denominator. This makes error values too large.
    #errorDT=np.sqrt(abs(sum(np.square(TimeShiftArray)*Correlations)/sum(Correlations)))

    
    #4)Visual CHECK
    if (plotSwitch!=0):
        maxCorrEstimate=np.arange(dt-errorDT,dt+errorDT)
        maxCorrValue=np.array([max(Correlations)]*len(maxCorrEstimate))
        plt.subplot(2,1,1)
        plt.plot(TimeShiftArray,Correlations,'g',linewidth=5)
        plt.plot(TimeShiftArray[zero],Correlations[zero],'b.',linewidth=5)
        plt.plot(maxCorrEstimate,maxCorrValue,'r',linewidth=6)
        plt.plot(theoryDT,max(Correlations),'k.',linewidth=8)
        plt.subplot(2,1,2)
        minT=pad
        maxT=len(Correlations)*TimeDelay
        ExperimentFrame=np.arange(minT,maxT)
        plt.plot(ExperimentFrame,bigA[minT:maxT],'b--',ExperimentFrame,bigB[minT:maxT],'r--')
        starterTime=np.arange(TimeStart,TimeStart+TimeBlock)
        expectedCorrTime=np.arange(TimeStart+dt,TimeStart+dt+TimeBlock)
        corrMiniB=bigB[TimeStart+dt:TimeStart+dt+TimeBlock]
        if (Debugger!=0):
            print 'starterTime: {}'.format(len(starterTime))
            print 'miniA: {}'.format(len(miniA))
            print 'expectedCorrTime: {}'.format(len(expectedCorrTime))
            print 'corrMiniB: {}'.format(len(corrMiniB))
        
        minis=plt.plot(starterTime,miniA,'g',expectedCorrTime,corrMiniB,'m')
        plt.setp(minis,linewidth=2.0)
    return dt,errorDT,percentCorrelation

def sampleArrays():
    A=0.0
    while A==0.0:
        A=(random.random()*100)-(random.random()*100)
    f=random.random()*100
    SNR=(random.random())*18+2
    primecheck=0
    while primecheck==0:
        length=int((random.random())*2000)+200
        primecheck=1
        for x in range(2,length-1):
            if ((length%x)==0):
                primecheck=0
    shift=int(random.random()*length/4)-int(random.random()*length/4)
    print 'Waveform is {}sin({}t).\n shift is {}.\n SNR= {}.\n Array length is {}'.format(A,f,shift,SNR,length)
    NoiseA=np.array([0.0]*length)
    NoiseB=np.array([0.0]*length)
    arrayA=np.array([0.0]*length)
    arrayB=np.array([0.0]*length)
    for x in range(0,length):
        NoiseA[x]=random.gauss(0,A/SNR)
        NoiseB[x]=random.gauss(0,A/SNR)
    for x in range(0,length):
        arrayA[x]=A*np.sin(f*x)
        arrayB[x]=A*np.sin(f*(x+shift))
    timeArray=np.arange(length)
    if (Debugger!=0):
        plt.plot(timeArray,arrayB,'r--',timeArray,NoiseB,'b--')
        plt.show()
    arrayA=arrayA+NoiseA
    arrayB=arrayB+NoiseB
    plt.subplot(2,1,1)
    arrays=plt.plot(timeArray,arrayA,'r--',timeArray,arrayB, 'b--')
    plt.setp(arrays,linewidth=2.0)
    return arrayA,arrayB

def WavePackets():
    print 'TDEBasic.wavepackets:'
    sineA,sineB,t,tShift,Amp,omega,shift=shiftedSines()
    length=len(sineA)
    print 'Base sinewave is sin({}t)\nshift is {}\nlength is {}'.format(omega,shift,length)
    packetA,packetB,sigma=Packets(t,tShift,length)
    print 'sigma of wavepacket is {}'.format(sigma)
    maxOfPacket=max([max(packetA),max(packetB)])
    Amp=Amp*maxOfPacket
    print 'max amplitude is {}'.format(Amp)
    noiseA,noiseB,SNR=NoiseMaker(length,Amp)
    print 'SNR is {}'.format(SNR)
    wavepacketA=(sineA*packetA)+noiseA
    wavepacketB=(sineB*packetB)+noiseB
    if(len(wavepacketA)!=length):
        print 'error in array math!'
    if (Debugger!=0):
        wavepacketPlotter(t,wavepacketA,wavepacketB)
        plt.show()
    print '\n'
    return wavepacketA,wavepacketB,shift
    

def shiftedSines():
    Amp=random.random()*100
    omega=2*pi*(random.random()*100.0)
    prime=0
    while prime==0:
        length=int(random.uniform(200,2000)+0.4999)
        for x in range(2,length-1):
            if (length%x)!=0:
                prime=1
    shift=int(random.random()*(length/4)+0.4999)-int(random.random()*(length/4)+0.4999)
    t=np.arange(length)
    t=t-int(0.499+(random.random()*length/2.0))-int((length/4.0)+.499)#Off shifting wave packet so the wave packet does not start at the very beginning of trace.
    tShift=t-shift
    sineA=Amp*np.sin(omega*t)
    sineB=Amp*np.sin(omega*tShift)
    return sineA,sineB,t,tShift,Amp,omega,shift

def NoiseMaker(length,Amp):
    SNR=(random.random()*18)+2
    #SNR=1000
    noiseA=[0.0]*length
    noiseB=[0.0]*length
    for x in range(0,len(noiseA)):
        noiseA[x]=random.gauss(0,Amp/SNR)
    for x in range(0,len(noiseB)):
        noiseB[x]=random.gauss(0,Amp/SNR)
    return noiseA,noiseB,SNR

def Packets(t,tShift,shift):
    packetLength=2*(random.random()*random.random()*shift/4)
    sigma=packetLength/2.0
    #normalization=1/(sigma*np.sqrt(2*pi))
    normalization=random.random()
    packetA=normalization*np.exp(-(t/packetLength)**2)
    packetB=normalization*np.exp(-(tShift/packetLength)**2)
    return packetA,packetB,sigma

def wavepacketPlotter(t,wavepacketA,wavepacketB):
    plt.subplot(3,1,1)
    arrays=plt.plot(t,wavepacketA,'b--',t,wavepacketB, 'r--')
    plt.setp(arrays,linewidth=2.0)
    tMaxA=np.argmax(wavepacketA)
    tMaxB=np.argmax(wavepacketB)
    plt.plot(t[tMaxA],max(wavepacketA),'go',t[tMaxB],max(wavepacketB),'go')

def Analyzer(shift,dt,Correlation,error):
    actualError=(dt-shift)/shift
    if ((dt-error<=shift) and (shift<=dt+error)):
        print "\n{}% error".format(actualError*100)
    else:
        print "\n{}% error".format(actualError*100)
        print "Actual shift is outside of errorbars"
    if Correlation<51:
        print "Low Correlation Value"

#For later
"""    
def multiPackets(dtx,dty):
    Amp=100*random.random()
    omega=2*pi*100*random.random()
    timelength=int(random.random()*1800)+200
    t=np.arange(timelength)
    t=t-int(timelength/4)-int(random.random()*timelength/2)
    sineA=Amp*np.sin(omega*t)
    packetA=true

def velocityScenario(dx,dy,vx,vy):
    shiftX=dx/vx
    shiftY=dy/vy
    return dtx,dty
"""
