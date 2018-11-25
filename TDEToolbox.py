"""
Compiliation of finalized tools used for TDE algorithms, including a more updated version of TDEBasic.TDECorrelate
"""
import numpy as np
import random
import matplotlib.pyplot as plt
import pylab
Debugger=0
pi=np.pi

def TDECorrelate(bigA,bigB,TimeStart,TimeBlock,StepSize,StepNumber):
    #1)Initialize
    bigA=np.array(bigA)-np.mean(bigA)
    bigB=np.array(bigB)-np.mean(bigB)
    TotalCorrCells=(2*StepNumber)+1
    Correlations=np.array([0.0]*TotalCorrCells)
    TimeShiftArray=np.array([0.0]*len(Correlations))
    pad=TimeStart-(StepSize*StepNumber)
    miniA=bigA[TimeStart:(TimeStart+TimeBlock)]
    miniA=miniA-np.mean(miniA)
    starter=pad
    
    #2) CORRELATIONS
    for x in range(0,len(Correlations)):
        TimeShiftArray[x]=starter-TimeStart
        if TimeShiftArray[x]==0:
            zero=x
        miniB=bigB[starter:(starter+TimeBlock)]
        miniB=miniB-np.mean(miniB)
        if (Debugger!=0):
            print(len(miniA)-len(miniB))
        Correlations[x]=np.sum(miniA*miniB)/(TimeBlock*(np.std(miniA))*(np.std(miniB)))
        del miniB
        starter=starter+StepSize
        
    #3)Analysis
    TimeShiftArray=np.array(TimeShiftArray)
    i=(np.argmax(Correlations))
    dt=TimeShiftArray[i]
    Correlation=max(Correlations)*100
    posiError=i+1
    negiError=i-1
    posiEndCounter=1
    negiEndCounter=1
    while (Correlations[posiError]>0.5*max(Correlations)):
        posiEndCounter=posiEndCounter+1
        posiError=posiError+1
        if (i==len(Correlations)):
            print "Boom! Crash! Dead."
    while (Correlations[negiError]>0.5*max(Correlations)):
        negiEndCounter=negiEndCounter+1
        negiError=negiError-1
        if (i==0):
            print "Boom! Crash! Dead."
    if (posiEndCounter>negiEndCounter):
        error=(posiEndCounter*StepSize)
    else:
        error=(negiEndCounter*StepSize)
    """
    if (min(Correlations)<0):
        NormCorre=Correlations+abs(min(Correlations))
    else:
        NormCorre=Correlations-abs(min(Correlations))
    NormCorreFactor=1/sum(NormCorre)
    NormCorre=NormCorreFactor*NormCorre
    varriance=sum((NormCorre*np.square(TimeShiftArray)))
    mean=sum((NormCorre*TimeShiftArray))
    error=np.sqrt(varriance-np.square(mean))

    normCorrelations=1/sum(np.square(Correlations))
    Correlations=normCorrelations*Correlations
    #error=np.sqrt(sum(np.square(Correlations*TimeShiftArray))/sum(np.square(Correlations)))
    error=np.sqrt(abs(sum(np.square(Correlations*TimeShiftArray))/len(Correlations)))
    """
    return dt,Correlation,error

def padder(length,TimeBlock,StepSize,StepNumber,TestPoints):
    tmin=StepSize*StepNumber
    tmax=length-(TimeBlock+(StepSize*StepNumber))
    timeStartList=np.linspace(tmin,tmax,TestPoints)
    for x in range(0,TestPoints):
        timeStartList[x]=int(timeStartList[x])
    if (Debugger!=0):
        print timeStartList
    return timeStartList

def dtFinder(arrayA,arrayB,TimeBlock,StepSize,StepNumber,TestPoints):
    timeStartList=padder(len(arrayA),TimeBlock,StepSize,StepNumber,TestPoints)
    maxCorrelation=0.0
    maxDT=0.0
    maxError=0.0
    for x in range(0,len(timeStartList)):
        TimeStart=timeStartList[x]
        dt,Correlation,error=TDECorrelate(arrayA,arrayB,TimeStart,TimeBlock,StepSize,StepNumber)
        if (Correlation>maxCorrelation):
            maxDT=dt
            maxCorrelation=Correlation
            maxError=error
    print "dt={}+/-{} {}% Correlation".format(maxDT,maxError,maxCorrelation)
    return maxDT,maxCorrelation,maxError

def VelocityMaker(dt,errorDT,x1,x2):
    V=(x2-x1)/dt
    errorV=abs(V*abs(errorDT/dt))
    print "V={}+/-{}".format(V,abs(errorV))
    return V,errorV

def printer(dt,error,correlation):
    print "{}+/-{} {}%Correlation".format(dt,error,correlation)

def dtAnalyzer(dt,error,correlation,correctDT):
    if (correctDT==0.0):
        correctDT=0.1
    percentError=100*(dt-correctDT)/correctDT
    if (correlation<=50):
        print "Low Correlation for dt"
    if (dt<(correctDT-error) or dt>(correctDT+error)):
        print "dt Calculation is outside of error bars"
    print "{}% error from correct dt Value".format(percentError)

def velocityAnalyzer(vMeasure,vTheory,vUncertainty):
    """
    print "vMeasure={}".format(vMeasure)
    print "vTheory={}".format(vTheory)
    print "vUncertainty={}".format(vUncertainty)
    """
    velocError=(vMeasure-vTheory)*100/vTheory
    print "{}% error from correct velocity value".format(velocError)
    if ((vMeasure<(vTheory-vUncertainty)) or (vMeasure>(vTheory+vUncertainty))):
        print "Velocity Calculation is outside of error bars"
    if (vUncertainty/vMeasure>1):
        print "Warning, velocity uncertainty is greater than 100%"
        
def Visualizer(bigA,bigB,measuredDT,t,TimeBlock,StepSize,StepNumber,first,theoryDT):
    #1)Initialize
    #print (np.argmin(abs(t)))
    TimeStart=(measuredDT+np.argmin(abs(t)))
    bigA=np.array(bigA)-np.mean(bigA)
    bigB=np.array(bigB)-np.mean(bigB)
    TotalCorrCells=(2*StepNumber)+1
    Correlations=np.array([0.0]*TotalCorrCells)
    TimeShiftArray=np.array([0.0]*len(Correlations))
    pad=TimeStart-(StepSize*StepNumber)
    miniA=bigA[TimeStart:(TimeStart+TimeBlock)]
    miniA=miniA-np.mean(miniA)
    starter=pad
    
    #2) CORRELATIONS
    for x in range(0,len(Correlations)):
        TimeShiftArray[x]=starter-TimeStart
        if TimeShiftArray[x]==0:
            zero=x
        miniB=bigB[starter:(starter+TimeBlock)]
        miniB=miniB-np.mean(miniB)
        if (Debugger!=0):
            print(len(miniA)-len(miniB))
        Correlations[x]=np.sum(miniA*miniB)/(TimeBlock*(np.std(miniA))*(np.std(miniB)))
        del miniB
        starter=starter+StepSize
        
    #3)Analysis
    TimeShiftArray=np.array(TimeShiftArray)
    i=(np.argmax(Correlations))
    dt=TimeShiftArray[i]
    Correlation=max(Correlations)*100
    """
    if (min(Correlations)<0):
        NormCorre=Correlations+abs(min(Correlations))
    else:
        NormCorre=Correlations-abs(min(Correlations))
    NormCorreFactor=1/sum(NormCorre)
    NormCorre=NormCorreFactor*NormCorre
    varriance=sum((NormCorre*np.square(TimeShiftArray)))
    mean=sum((NormCorre*TimeShiftArray))
    error=np.sqrt(varriance-np.square(mean))
    """
    posiError=i+1
    negiError=i-1
    posiEndCounter=1
    negiEndCounter=1
    while (Correlations[posiError]>0.5*max(Correlations)):
        posiEndCounter=posiEndCounter+1
        posiError=posiError+1
    while (Correlations[negiError]>0.5*max(Correlations)):
        negiEndCounter=negiEndCounter+1
        negiError=negiError-1
    if (posiEndCounter>negiEndCounter):
        error=(posiEndCounter*StepSize)
    else:
        error=(negiEndCounter*StepSize)
    """
    if (min(Correlations)<0):
        NormCorre=Correlations+abs(min(Correlations))
    else:
        NormCorre=Correlations-abs(min(Correlations))
    NormCorreFactor=1/sum(NormCorre)
    NormCorre=NormCorreFactor*NormCorre
    error=np.sqrt(abs(sum(NormCorre*np.square(TimeShiftArray))))
    """
    #4)Visual Check
    maxCorrEstimate=np.arange(dt-error,dt+error)
    maxCorrValue=np.array([max(Correlations)]*len(maxCorrEstimate))
    plt.subplot(4,1,first)
    plt.plot(TimeShiftArray,Correlations,'g',linewidth=5)
    plt.plot(TimeShiftArray[zero],Correlations[zero],'b.',linewidth=5)
    plt.plot(maxCorrEstimate,maxCorrValue,'r',linewidth=6)
    plt.plot(theoryDT,max(Correlations),'k.',linewidth=8)
    second=first+1
    plt.subplot(4,1,second)
    minT=pad
    maxT=len(Correlations)*StepSize+pad
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

def MaxChecker(wave1,wave2,x1,x2):
    dt=np.argmax(wave2)-np.argmax(wave1)
    if (dt==0):
        V=0
    else:
        V=(x2-x1)/dt
    print "MaxChecker states dt={}".format(dt)
    print "MaxChecker states V={}".format(V)
