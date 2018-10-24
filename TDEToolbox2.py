"""
A new version TDEToolbox. Including Correlation measurement seperated from error measurement, a class generator making pixel objects
"""
import numpy as np
import random
import matplotlib.pyplot as plt
import pylab
Debugger=0
pi=np.pi

def Measurer(pixelA,pixelB,numberOfRuns,Visualizer,blockSize,stepSize,stepNumber,testPoints):
    x1=pixelA.xCoor
    y1=pixelA.yCoor
    x2=pixelB.xCoor
    y2=pixelB.yCoor
    
    for N in range(0,numberOfRuns):
        aN=pixelA.timeData[N,:]
        bN=pixelB.timeData[N,:]
        if (Visualizer!=0):
            dt,Correlation,CorrelationsN,tStart=dtFinder(aN,bN,blockSize,stepSize,stepNumber,testPoints,Visualizer)
        else:
            dt,Correlation=dtFinder(aN,bN,blockSize,stepSize,stepNumber,testPoints,Visualizer)
        if (N==0):
            dtB=[dt]
            CorrelationB=[Correlation]
            if (Visualizer!=0):
                Correlations=CorrelationsN
                tStartB=[tStart] #Visualizer Only
        else:
            dtB.append(dt)
            CorrelationB.append(Correlation)
            if (Visualizer!=0):
                Correlations=np.vstack((Correlations,CorrelationsN))
                tStartB.append(tStart) #Visualizer Only
    
    dtB,errorDTB=meanAndDeviation(dtB)
    errorDTB=errorDTB+stepSize/2.0 #Note: added in because sqrt(<(dt)^2>)==0 when calculating wavepackets. 
    CorrelationB,errorCorrelationB=meanAndDeviation(CorrelationB)
    Vx,Vy,errorVX,errorVY=velocityMaker(x1,y2,x2,y2,dtB,errorDTB)
    if (Visualizer!=0):
        tStart=int(np.mean(tStartB)+0.499) #Visualizer Only
    
    pixelB.dtAndCorrelation(dtB,errorDTB,CorrelationB,errorCorrelationB)
    pixelB.velocityRecorder(Vx,Vy,errorVX,errorVY)
    #pixelB.Printer() #<--In real case, remove
    if (Visualizer!=0):
        aveCorrelations=np.mean(Correlations,axis=0)#Visualizer only
        stdCorrelations=np.std(Correlations,axis=0)#Visualizer only
    if (Visualizer!=0):
        return pixelB,aveCorrelations,stdCorrelations,tStart
    else:
        return pixelB

def TDECorrelate(bigA,bigB,TimeStart,TimeBlock,StepSize,StepNumber,Visualizer):
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
    #4)Visualizer Switch
    if (Visualizer==0):
        return dt,Correlation #For final version
    else:
        return dt,Correlation,Correlations #For testing purposes

def dtFinder(arrayA,arrayB,TimeBlock,StepSize,StepNumber,TestPoints,Visualizer):
    timeStartList=padder(len(arrayA),TimeBlock,StepSize,StepNumber,TestPoints)
    maxCorrelation=0.0
    maxCorrelations=[0.0]*(2*StepNumber+1)
    maxDT=0.0
    for x in range(0,len(timeStartList)):
        TimeStart=timeStartList[x]
        if (Visualizer==0):
            dt,Correlation=TDECorrelate(arrayA,arrayB,TimeStart,TimeBlock,StepSize,StepNumber,Visualizer)
        else:
            dt,Correlation,Correlations=TDECorrelate(arrayA,arrayB,TimeStart,TimeBlock,StepSize,StepNumber,Visualizer)
        if (Correlation>maxCorrelation):
            maxDT=dt
            maxCorrelation=Correlation
            if (Visualizer!=0):
                maxCorrelations=Correlations
                tStart=TimeStart
    if (Visualizer!=0):
        return maxDT,maxCorrelation,maxCorrelations,tStart
    else:
        return maxDT,maxCorrelation

def padder(length,TimeBlock,StepSize,StepNumber,TestPoints):
    tmin=StepSize*StepNumber
    tmax=length-(TimeBlock+(StepSize*StepNumber))
    timeStartList=np.linspace(tmin,tmax,TestPoints)
    for x in range(0,TestPoints):
        timeStartList[x]=int(timeStartList[x])
    if (Debugger!=0):
        print timeStartList
    return timeStartList

def meanAndDeviation(a):
    return np.mean(a),np.std(a)

def velocityMaker(x1,y1,x2,y2,dt,errorDT):
    dx=x2-x1
    dy=y2-y1
    if (dt==0):
        dt=0.001
    if (dx==0):
        Vx=0.000
    else:
        Vx=dx/dt
    if (dy==0):
        Vy=0.000
    else:
        print dy
        Vy=dy/dt
    errorVX=abs(Vx*(errorDT/dt))
    errorVY=abs(Vy*(errorDT/dt))
    return Vx,Vy,errorVX,errorVY

"""(Re-written in TDETestCases due object oriented code)
def dtAnalyzer(Correlation,errorCorre,measuredDT,errorDT,theoryDT,maxDT):
    if (Correlation<51):
        print "Low Correlation"
    if (errorCorre>=Correlation):
        print "Imprecise Correlation"
    print "dt measurement is {}% from theoryDT".format((measuredDT-theoryDT)*(100.0/theoryDT))
    if (theoryDT<(measuredDT-errorDT) or theoryDT>(measuredDT+errorDT)):
        print "theoryDT is outside of errorbars"
    print "dt measurement is {}% from maxDT".format((measuredDT-maxDT)*(100.0/maxDT))
    if (maxDT<(measuredDT-errorDT) or maxDT>(measuredDT+errorDT)):
        print "maxDT is outside of errorbars"
"""

"""(Re-written in TDETestCases due object oriented code)
def velocityAnalyzer(measuredVX,errorVX,theoryVX,maxVX,measuredVY,errorVY,theoryVY,maxVY):
    print "Vx measurement is {}% from theoryVX".format((measuredVX-theoryVX)*(100.0/theoryVX))
    if (theoryVX<(measuredVX-errorVX) or theoryVX>(measuredVX+errorVX)):
        print "theoryVX is outside of errorbars"
    print "Vx measurement is {}% from maxVX".format((measuredVX-maxVX)*(100.0/maxVX))
    if (maxVX<(measuredVX-errorVX) or maxVX>(measuredVX+errorVX)):
        print "maxVX is outside of errorbars"
    
    print "Vy measurement is {}% from theoryVY".format((measuredVY-theoryVY)*(100.0/theoryVY))
    if (theoryVY<(measuredVY-errorVY) or theoryVY>(measuredVY+errorVY)):
        print "theoryVY is outside of errorbars"
    print "Vy measurement is {}% from maxVY".format((measuredVY-maxVY)*(100.0/maxVY))
    if (maxVY<(measuredVY-errorVY) or maxVY>(measuredVY+errorVY)):
        print "maxVY is outside of errorbars"
"""
    
def MaxChecker(wave1,wave2,x1,y1,x2,y2):
    dt=np.argmax(wave2)-np.argmax(wave1)
    dx=x2-x1
    dy=y2-y1
    if (dx==0.0 or dt==0.0):
        Vx=0.0
    else:
        Vx=dx/dt
    if (dy==0.0 or dt==0.0):
        Vy=0.0
    else:
        Vy=dy/dt
    return dt,Vx,Vy

def MaxChecker2(matrix1,matrix2,x1,y1,x2,y2):
    #Like MaxChecker,but mainly used for repeatable tests.
    #Will give expected value(mean(maxes)) and varied value(std((max))
    if (Debugger!=0):
        print sum(sum(matrix1-matrix2))
    repeatables=len(matrix1[:,0])
    dtB=[0.0]*repeatables
    VxB=[0.0]*repeatables
    VyB=[0.0]*repeatables
    
    for cell in range(0,repeatables):
        wave1=matrix1[cell,:]
        wave2=matrix2[cell,:]
        maxDT,maxVX,maxVY=MaxChecker(wave1,wave2,x1,y1,x2,y2)
        dtB[cell]=maxDT
        VxB[cell]=maxVX
        VyB[cell]=maxVY
        del wave1
        del wave2
    dt,errorDT=meanAndDeviation(dtB)
    Vx,errorVX=meanAndDeviation(VxB)
    Vy,errorVY=meanAndDeviation(VyB)
    return dt,Vx,Vy,errorDT,errorVX,errorVY
    

def Analyzer(pointA,pointB,theoryDT,theoryVX,theoryVY):
    matrix1=pointA.timeData
    matrix2=pointB.timeData
    x1=pointA.xCoor
    y1=pointA.yCoor
    x2=pointB.xCoor
    y2=pointB.yCoor
    
    maxDT,maxVX,maxVY,maxErrorDT,maxErrorVX,maxErrorVY=MaxChecker2(matrix1,matrix2,x1,y1,x2,y2)
    #dt=pointB.dt
    #errorDT=pointB.errorDT
    print "\nfor ({},{})".format(x2,y2)
    print "Theory dt={}. maxChecker states dt={}+/-{}".format(theoryDT,maxDT,maxErrorDT)
    print "Theory Vx={}. maxChecker states Vx={}+/-{}".format(theoryVX,maxVX,maxErrorVX)
    print "Theory Vy={}. maxChecker states Vy={}+/-{}".format(theoryVY,maxVY,maxErrorVY)
    pointB.dtAnalyzer(theoryDT,maxDT,maxErrorDT)
    pointB.velocityAnalyzer(theoryVX,maxVX,maxErrorVX,theoryVY,maxVY,maxErrorVY)
    
def Visualizer(pointA,pointB,aveC,aveEC,tStart,theoryDT,blockSize,stepSize,stepNumber,first):
    reach=stepSize*stepNumber
    second=first+1
    matrix1=pointA.timeData
    matrix2=pointB.timeData
    x1=pointA.xCoor
    y1=pointA.yCoor
    x2=pointB.xCoor
    y2=pointB.yCoor
    wave1=pointA.averageData
    wave2=pointB.averageData
    maxDT,maxVX,maxVY,maxEDT,maxEVX,maxEVY=MaxChecker2(matrix1,matrix2,x1,y1,x2,y2)
    timeShiftArray=(np.arange(2*stepNumber+1)*stepSize)-reach
    eRange=np.arange(pointB.dt-pointB.errorDT,pointB.dt+pointB.errorDT)
    ExperimentFrame=np.arange(tStart-reach,tStart+reach+stepSize)
    ZeroFrame=np.arange(tStart,tStart+blockSize)
    maxErrorFrame=np.arange((maxDT-maxEDT),(maxDT+maxEDT))
    dtFrame=np.array(ZeroFrame)+int(pointB.dt+0.499)
    cell=0
    zero=-1
    while (zero<0):
        if (timeShiftArray[cell]==0):
            zero=cell
        else:
            cell=cell+1
    
    plt.subplot(4,1,first)
    #plt.errorbar(timeShiftArray,aveC,aveEC,'g',linewidth=5)
    plt.errorbar(timeShiftArray,aveC,aveEC,color='g',ecolor='g')
    plt.plot(timeShiftArray[zero],aveC[zero],'b.',linewidth=5)
    plt.plot(pointB.dt,max(aveC),'rs',linewidth=5)
    plt.plot(eRange,[max(aveC)]*len(eRange),'rs',linewidth=5)
    plt.plot(theoryDT,max(aveC),'ko',linewidth=8)
    plt.plot(maxErrorFrame,[max(aveC)]*len(maxErrorFrame),'yo',linewidth=10)

    plt.subplot(4,1,second)
    plt.plot(ExperimentFrame,wave1[ExperimentFrame],'b--',ExperimentFrame,wave2[ExperimentFrame],'r--')
    miniA=wave1[tStart:(tStart+blockSize)]
    miniB=wave2[dtFrame]
    minis=plt.plot(ZeroFrame,miniA,'g',dtFrame,miniB,'m')
    plt.setp(minis,linewidth=2.0)

def subplotFinishier(xAxis,yAxis,title):
    plt.xlabel(xAxis,fontsize=15)
    plt.ylabel(yAxis,fontsize=15)
    plt.title(title,fontsize=20)
    
def StepSizeVisual(subFig,stepSize,velocity,errorVelocity,maxVelocity):
    plt.subplot(4,1,subFig)
    measurements=plt.errorbar(stepSize,velocity,errorVelocity,color='g',ecolor='g',linewidth=2,label="measurements")
    expectation=plt.plot(stepSize,[maxVelocity]*len(stepSize),'y',linewidth=2,label="expectation")
    if (subFig==1):
        subtitle="Vx for X Case"
    elif (subFig==2):
        subtitle="Vy for X Case"
    elif (subFig==3):
        subtitle="Vx for Y Case"
    elif (subFig==4):
        subtitle="Vy for Y Case"
    else:
        print "error"
    plt.legend(loc='upper right')
    subplotFinishier("stepsize","velocity",subtitle)
    
def StepSizeVisualDT(subFig,stepSize,dt,errorDT,maxDT,maxErrorDT):
    plt.subplot(2,1,subFig)
    measurements=plt.errorbar(stepSize,dt,errorDT,color='g',ecolor='g',linewidth=2,label="measurements")
    expectation=plt.errorbar(stepSize,[maxDT]*len(stepSize),[maxErrorDT]*len(stepSize),color='y',ecolor='y',linewidth=2,label="expectation")
    if (subFig==1):
        subtitle="dt for X Case"
    elif (subFig==2):
        subtitle="dt for Y Case"
    plt.legend(loc='upper right')
    subplotFinishier("stepsize","dt",subtitle)
