"""
TDETest
Let's Test the TDE calculations
"""
import TDEBasic as TDE
import TDEToolbox as TBox
import TDETestCases as TCases
import matplotlib.pyplot as plt
import TDEToolbox2 as TBox2
#import random
import numpy as np
Debugger=0
TestNumber=5

#TEST 1: TDE CALCULATION W/ USER INPUTED TimeStart <Successful>
if (TestNumber==1):
    arrayA,arrayB,shift=TDE.WavePackets()
    if (shift==0):
        shift=0.001
    print 'Maximum of first waveform at t={}'.format(np.argmax(arrayA))
    TimeStart=input("Input desired starting time: ")
    TimeBlock=input("Input desired timeblock: ")
    TimeDelay=input("Input desired delay shift: ")
    dt,error,Correlation=TDE.TDECorrelate(arrayA,arrayB,TimeStart,TimeBlock,TimeDelay,1,shift)
    print "dt={}+/-{}. {}% correlation".format(dt,error,Correlation)
    TDE.Analyzer(shift,dt,Correlation,error)
    plt.show()

#Test 2: TDE CALCULATION W/O USER INPUTED TimeStart <good>
if (TestNumber==2):
    arrayA,arrayB,shift=TDE.WavePackets()
    if (shift==0):
        shift=0.001
    TimeBlock=input("Input desired timeblock: ")
    StepSize=input("Input desired stepSize: ")
    trials=int((len(arrayA)-TimeBlock)/StepSize)
    dtEstimates=[0.0]*trials
    errorEstimates=[0.0]*trials
    correlationEstimates=[0.0]*trials
    for x in range(0,trials):
        TimeStart=x*StepSize
        dt,errorDT,Correlation=TDE.TDECorrelate(arrayA,arrayB,TimeStart,TimeBlock,StepSize,0,shift)
        correlationEstimates[x]=Correlation
    i=np.argmax(correlationEstimates)
    dt,errorDT,Correlation=TDE.TDECorrelate(arrayA,arrayB,(i*StepSize),TimeBlock,StepSize,1,shift)

    print "dt={}+/-{}. {}% correlation".format(dt,errorDT,Correlation)
    TDE.Analyzer(shift,dt,Correlation,errorDT)
    plt.show()

#Test 3) 3 point generator <Failure>
#Abandoned due to misconceptions of error analysis. Running this test will lead to errors
if (TestNumber==3):
    origin,dxPoint,dyPoint,Vx,Vy,x,y,t=TCases.ThreePointGenerator()
    if (Debugger!=0):
        print t
    expectedDTx=x/Vx
    print "expectedDTx={}".format(expectedDTx)
    expectedDTy=y/Vy
    print "expectedDTy={}".format(expectedDTy)
    length=len(origin)
    TimeBlock=200
    reach=250
    StepSize=5
    StepNumber=int(reach/StepSize)
    TestPoints=50
    
    print "\nX test"
    TBox.MaxChecker(origin,dxPoint,0,x)
    print "Theory dtX={}".format(expectedDTx)
    print "Theory Vx={}".format(Vx)
    dtX,CoorX,ErrorX=TBox.dtFinder(origin,dxPoint,TimeBlock,StepSize,StepNumber,TestPoints)
    if (dtX==0.0):
        print "dtX=0"
        dtX=1
    measuredVx,measuredErrorVx=TBox.VelocityMaker(dtX,ErrorX,0,x)
    TBox.dtAnalyzer(dtX,ErrorX,CoorX,expectedDTx)
    TBox.velocityAnalyzer(measuredVx,Vx,measuredErrorVx)

    print "\nY test"
    TBox.MaxChecker(origin,dyPoint,0,y)
    print "Theory dtY={}".format(expectedDTy)
    print "Theory Vy={}".format(Vy)
    dtY,CoorY,ErrorY=TBox.dtFinder(origin,dyPoint,TimeBlock,StepSize,StepNumber,TestPoints)
    if (dtY==0.0):
        print "dtY=0"
        dtY=1
    measuredVy,measuredErrorVy=TBox.VelocityMaker(dtY,ErrorY,0,y)
    TBox.dtAnalyzer(dtY,ErrorY,CoorY,expectedDTy)
    TBox.velocityAnalyzer(measuredVy,Vy,measuredErrorVy)
    
    TBox.Visualizer(origin,dxPoint,int(dtX),t,TimeBlock,StepSize,StepNumber,1,expectedDTx)
    TBox.Visualizer(origin,dyPoint,int(dtY),t,TimeBlock,StepSize,StepNumber,3,expectedDTy)
    plt.show()
        

#Test 4) 3 point generator <inconclusive>
if (TestNumber==4):
    Axis=2
    #arrayA,arrayB,shift=TDE.WavePackets()
    origin,dxPoint,dyPoint,Vx,Vy,x,y,t=TCases.ThreePointGenerator()
    expectedDTx=x/Vx
    print "expectedDTx={}".format(expectedDTx)
    expectedDTy=y/Vy
    print "expectedDTy={}".format(expectedDTy)
    length=len(origin)

    arrayA=origin
    if (Axis==1):
        print "X value test"
        arrayB=dxPoint
        shift=expectedDTx
    elif (Axis==2):
        print "Y value test"
        arrayB=dyPoint
        shift=expectedDTy
    if (shift==0):
        shift=0.001
    TimeBlock=input("Input desired timeblock: ")
    StepSize=input("Input desired stepSize: ")
    trials=int((len(arrayA)-TimeBlock)/StepSize)
    dtEstimates=[0.0]*trials
    errorEstimates=[0.0]*trials
    correlationEstimates=[0.0]*trials
    for x in range(0,trials):
        TimeStart=x*StepSize
        dt,errorDT,Correlation=TDE.TDECorrelate(arrayA,arrayB,TimeStart,TimeBlock,StepSize,0,shift)
        correlationEstimates[x]=Correlation
    i=np.argmax(correlationEstimates)
    dt,errorDT,Correlation=TDE.TDECorrelate(arrayA,arrayB,(i*StepSize),TimeBlock,StepSize,1,shift)

    print "dt={}+/-{}. {}% correlation".format(dt,errorDT,Correlation)
    TDE.Analyzer(shift,dt,Correlation,errorDT)
    plt.show()

#Test 5: Creating an error bar via repetitions <Satisfactory>
#If procedure ran in real life, it would need N number of runs
#NB: Measuring algorithm is here, not in TBox2.py (TBox2.Measurer is only for Tests 6 & 7)
if (TestNumber==5):
    #Input Stage (User inputted or user scripted)
    blockSize=200
    reach=200
    stepSize=5
    stepNumber=int(reach/stepSize)
    testPoints=25 #Deterimines number of timeStarts
    numberOfRuns=20
    Visualizer=1 # 1=on 0=off
    
    #Data
    originData,xData,yData,Vx0,Vy0,x,y,t=TCases.ThreePoint_NRunGenerator(numberOfRuns)
    theoryDTX0=x/Vx0
    theoryDTY0=y/Vy0
    length=len(originData[0,:])

    #Storage for each point
    point0=TCases.pixel(0,0,originData)
    point1=TCases.pixel(x,0,xData)
    point2=TCases.pixel(0,y,yData)
        
    #Measurement Stage (Goal: get this generalized as possible for future use)
    x0=point0.xCoor
    y0=point0.yCoor
    for point in range(1,3):
        if (point==1):
            bData=point1.timeData
            xCoor=point1.xCoor
            yCoor=point1.yCoor
        elif (point==2):
            bData=point2.timeData
            xCoor=point2.xCoor
            yCoor=point2.yCoor
        else:
            print "error"
            break
        for N in range(0,numberOfRuns):
            aN=point0.timeData[N,:]
            bN=bData[N,:]
            if (Visualizer!=0):
                dt,Correlation,CorrelationsN,tStart=TBox2.dtFinder(aN,bN,blockSize,stepSize,stepNumber,testPoints,Visualizer)
            else:
                dt,Correlation=TBox2.dtFinder(aN,bN,blockSize,stepSize,stepNumber,testPoints,Visualizer)
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
        dtB,errorDTB=TBox2.meanAndDeviation(dtB)
        errorDTB=errorDTB+stepSize/2.0 #Note: added in because sqrt(<(dt)^2>)==0 when calculating wavepackets.
        CorrelationB,errorCorrelationB=TBox2.meanAndDeviation(CorrelationB)
        Vx,Vy,errorVX,errorVY=TBox2.velocityMaker(x0,y0,xCoor,yCoor,dtB,errorDTB)
        if (Visualizer!=0):
            tStart=int(np.mean(tStartB)+0.499) #Visualizer Only
        
        if (point==1):
            point1.dtAndCorrelation(dtB,errorDTB,CorrelationB,errorCorrelationB)
            point1.velocityRecorder(Vx,Vy,errorVX,errorVY)
            point1.Printer() #<--In real case, remove
            if (Visualizer!=0):
                CorrelationsX=Correlations #Visualizer only
                tStartX=tStart #Visualizer only
        elif(point==2):
            point2.dtAndCorrelation(dtB,errorDTB,Correlation,errorCorrelationB)
            point2.velocityRecorder(Vx,Vy,errorVX,errorVY)
            point2.Printer() #<--In real case, remove
            if (Visualizer!=0):
                CorrelationsY=Correlations #Visualizer only
                tStartY=tStart #Visualizer only
        else:
            print "error"
            break
        del bData
        
    #Analysis Stage (Only for testing purposes)
    wave0=point0.averageData
    matrix0=point0.timeData
    x1=point0.xCoor
    y1=point0.yCoor
    for Test in range(1,3):
        if (Test==1):
            theoryDT=theoryDTX0
            theoryVX=Vx0
            theoryVY=0.0
            TBox2.Analyzer(point0,point1,theoryDT,theoryVX,theoryVY)
        elif (Test==2):
            theoryDT=theoryDTY0
            theoryVX=0.0
            theoryVY=Vy0
            TBox2.Analyzer(point0,point2,theoryDT,theoryVX,theoryVY)
        else:
            print "error"
            break
        
    #Visualizer (Only for testing purposes)
    aveCorrelationsX=np.mean(CorrelationsX,axis=0)
    stdCorrelationsX=np.std(CorrelationsX,axis=0)
    aveCorrelationsY=np.mean(CorrelationsY,axis=0)
    stdCorrelationsY=np.std(CorrelationsY,axis=0)
    TBox2.Visualizer(point0,point1,aveCorrelationsX,stdCorrelationsX,tStartX,theoryDTX0,blockSize,stepSize,stepNumber,1)
    TBox2.Visualizer(point0,point2,aveCorrelationsY,stdCorrelationsY,tStartY,theoryDTY0,blockSize,stepSize,stepNumber,3)
    plt.show()
   
#Test 6) StepSize Alterations <Inconclusive>
"""
Pretty similar to Test 5 except:
-Measurement Stage is now a function in TDEToolbox2.py
-Visualizer and Analyzer are turned off
-Each measurement case will be compared to maxChecker
"""

if (TestNumber==6):
    #Unchangeable Inputs
    blockSize=200
    reach=200
    testPoints=25 #Deterimines number of timeStarts
    numberOfRuns=5
    smallestStep=5
    LargestStep=20
    Visualizer=0 # 1=on 0=off
    
    #Data
    originData,xData,yData,Vx0,Vy0,x,y,t=TCases.ThreePoint_NRunGenerator(numberOfRuns)
    theoryDTX0=x/Vx0
    theoryDTY0=y/Vy0
    length=len(originData[0,:])

    #Storage for each point
    point0=TCases.pixel(0,0,originData)
    pointX=TCases.pixel(x,0,xData)
    pointY=TCases.pixel(0,y,yData)
    wave0=point0.averageData
    wave1=pointX.averageData
    wave2=pointY.averageData

    
    stepSizes=np.linspace(smallestStep,LargestStep,10)
    for size in range(0,len(stepSizes)):
        stepSizes[size]=int(stepSizes[size]+0.499)
    
    VxXCase=[0.0]*len(stepSizes) 
    VyXCase=[0.0]*len(stepSizes)
    VxYCase=[0.0]*len(stepSizes)
    VyYCase=[0.0]*len(stepSizes)
    
    errorVxXCase=[0.0]*len(stepSizes)
    errorVyXCase=[0.0]*len(stepSizes)
    errorVxYCase=[0.0]*len(stepSizes)
    errorVyYCase=[0.0]*len(stepSizes)
    
    dt1,Vx1,Vy1=TBox2.MaxChecker(wave0,wave1,point0.xCoor,point0.yCoor,pointX.xCoor,pointX.yCoor)
    dt2,Vx2,Vy2=TBox2.MaxChecker(wave0,wave2,point0.xCoor,point0.yCoor,pointY.xCoor,pointY.yCoor)
    maxVelocity=[Vx1,Vy1,Vx2,Vy2]
    
    i=0
    for stepSize in stepSizes:
        stepNumber=int(reach/stepSize)
        temporaryPoint1=TBox2.Measurer(point0,pointX,numberOfRuns,Visualizer,blockSize,stepSize,stepNumber,testPoints)
        temporaryPoint2=TBox2.Measurer(point0,pointY,numberOfRuns,Visualizer,blockSize,stepSize,stepNumber,testPoints)
        
        VxXCase[i]=temporaryPoint1.Vx
        VyXCase[i]=temporaryPoint1.Vy
        VxYCase[i]=temporaryPoint2.Vx
        VyYCase[i]=temporaryPoint2.Vy

        errorVxXCase[i]=temporaryPoint1.errorVX
        errorVyXCase[i]=temporaryPoint1.errorVY
        errorVxYCase[i]=temporaryPoint2.errorVX
        errorVyYCase[i]=temporaryPoint2.errorVY
        i+=1
        
    Velocities=np.vstack((VxXCase,VyXCase,VxYCase,VyYCase))
    errorVelocities=np.vstack((errorVxXCase,errorVyXCase,errorVxYCase,errorVyYCase))
    #print "stepSizes[{}],VelocitiesN[{}],errorVelocities[{}]".format(len(stepSizes),len(Velocities[0,:]),len(errorVelocities[0,:]))
    for subFig in range(1,5):
        TBox2.StepSizeVisual(subFig,stepSizes,Velocities[subFig-1,:],errorVelocities[subFig-1,:],maxVelocity[subFig-1])
    plt.show()

#Test 7 StepSize Alterations part 2: checking dt
"""
Like test 6, but examining the dt instead of the velocities.
"""
if (TestNumber==7):
    #unchangable inputs
    blockSize=200
    reach=400
    testPoints=25 #Deterimines number of timeStarts
    numberOfRuns=25
    smallestStep=5
    LargestStep=20
    Visualizer=0 # 1=on 0=off
    
    #Data
    originData,xData,yData,Vx0,Vy0,x,y,t=TCases.ThreePoint_NRunGenerator(numberOfRuns)
    theoryDTX0=x/Vx0
    theoryDTY0=y/Vy0
    length=len(originData[0,:])

    #Storage for each point
    point0=TCases.pixel(0,0,originData)
    pointX=TCases.pixel(x,0,xData)
    pointY=TCases.pixel(0,y,yData)
    wave0=point0.averageData
    wave1=pointX.averageData
    wave2=pointY.averageData
    matrix0=point0.timeData
    matrix1=pointX.timeData
    matrix2=pointY.timeData

    
    stepSizes=np.linspace(smallestStep,LargestStep,10)
    for size in range(0,len(stepSizes)):
        stepSizes[size]=int(stepSizes[size]+0.499)

    dtXCase=[0.0]*len(stepSizes)
    errorDTXCase=[0.0]*len(stepSizes)
    dtYCase=[0.0]*len(stepSizes)
    errorDTYCase=[0.0]*len(stepSizes)

    dt1,Vx1,Vy1,errorDT1,errorVX1,errorVY1=TBox2.MaxChecker2(matrix0,matrix1,point0.xCoor,point0.yCoor,pointX.xCoor,pointX.yCoor)
    dt2,Vx2,Vy2,errorDT2,errorVX2,errorVY2=TBox2.MaxChecker2(matrix0,matrix2,point0.xCoor,point0.yCoor,pointY.xCoor,pointY.yCoor)
    maxDTs=[dt1,dt2]
    maxErrorDTs=[errorDT1,errorDT2]
    i=0
    for stepSize in stepSizes:
        stepNumber=int(reach/stepSize)
        temporaryPoint1=TBox2.Measurer(point0,pointX,numberOfRuns,Visualizer,blockSize,stepSize,stepNumber,testPoints)
        temporaryPoint2=TBox2.Measurer(point0,pointY,numberOfRuns,Visualizer,blockSize,stepSize,stepNumber,testPoints)
        
        dtXCase[i]=temporaryPoint1.dt
        dtYCase[i]=temporaryPoint2.dt

        errorDTXCase[i]=temporaryPoint1.errorDT
        errorDTYCase[i]=temporaryPoint2.errorDT
        i=i+1
        del temporaryPoint1
        del temporaryPoint2
    
    dts=np.vstack((dtXCase,dtYCase))
    errorDTs=np.vstack((errorDTXCase,errorDTYCase))
    for subFig in range(1,3):
        TBox2.StepSizeVisualDT(subFig,stepSizes,dts[subFig-1,:],errorDTs[subFig-1,:],maxDTs[subFig-1],maxErrorDTs[subFig-1])
    plt.show()
    

    """
        if (Visualizer!=0):
            #Visualizer (Only for testing purposes)
            aveCorrelationsX=np.mean(CorrelationsX,axis=0)
            stdCorrelationsX=np.std(CorrelationsX,axis=0)
            aveCorrelationsY=np.mean(CorrelationsY,axis=0)
            stdCorrelationsY=np.std(CorrelationsY,axis=0)
            timeShiftArray=(np.arange(2*stepNumber+1)*stepSize)-reach
            ExperimentFrame=[tStart-reach:tStart+reach+stepSize]
            ZeroFrame=[tStart:tStart+blockSize]
            dtFrame=np.array(ZeroFrame)+dt
            cell=0
            zero=-1
            while zero<0:
                if (timeShiftArray[cell]==0):
                    zero=cell
            for sPlot in range(1,3):
                if (sPlot==1):
                    first=1
                    aveC=aveCorrelationsX
                    aveEC=stdCorrelationsX
                    eRange=errorRange1
                    theoryT=theory1
                    maxT=max1
                elif (sPlot==2):
                    first=3
                    aveC=aveCorrelationsY
                    aveEC=stdCorrelationsY
                    eRange=errorRange2
                    theoryT=theory2
                    maxT=max2
            second=first+1
        
            plt.subplot(4,1,first)
            plt.errorbar(TimeShiftArray,aveC,aveEC,'g',linewidth=5)
            plt.plot(TimeShiftArray[zero],aveC[zero],'b.',linewidth=5)
            plt.plot(eRange,[max(aveC)]*len(eRange),'r',linewidth=6)
            plt.plot(theoryT,max(aveC),'k.',linewidth=8)
            plt.plot(maxT,max(aveC),'y',linewidth=8)
    
            plt.subplot(4,1,second)
            plt.plot(ExperimentFrame,bigA[:maxT],'b--',ExperimentFrame,bigB[minT:maxT],'r--')
            starterTime=np.arange(TimeStart,TimeStart+TimeBlock)
            expectedCorrTime=np.arange(TimeStart+dt,TimeStart+dt+TimeBlock)
            corrMiniB=bigB[TimeStart+dt:TimeStart+dt+TimeBlock]
    """
    """
    originWave,dxWave,dyWave,length,xShift,yShift=TCases.ThreePointGenerator()
    TimeBlock=100
    StepSize=10
    StepNumber=30
    TestPoints=50
    print "X test"
    print "xShift={}".format(xShift)
    dtX,correlationX,errorX=TBox.dtFinder(originWave,dxWave,TimeBlock,StepSize,StepNumber,TestPoints)
    TBox.printer(dtX,errorX,correlationX)
    TBox.analyzer(dtX,errorX,correlationX,xShift)

    print "Y test"
    print "yShift={}".format(yShift)
    dtY,correlationY,errorY=TBox.dtFinder(originWave,dyWave,TimeBlock,StepSize,StepNumber,TestPoints)
    TBox.printer(dtY,errorY,correlationY)
    TBox.analyzer(dtY,errorY,correlationY,yShift)
    TBox.Visualizer(originWave,dxWave,dtX,TimeBlock,StepSize,StepNumber,1)
    TBox.Visualizer(originWave,dyWave,dtY,TimeBlock,StepSize,StepNumber,3)
    plt.show()
    """
"""
#Test ??? (Abandoned)  Multi-timepoint TDE Calculation (Coded Inputs)
TimeBlock=100
StepSize=10
StepNumber=5
Length=1999 #Because its the largest prime number under  2000
TestTimes=TBox.padder(Length,TimeBlock,StepSize,StepNumber)
arrayA,arrayB=TCases.___(length)
dt=[0.0]*len(TestTimes)
errors=[0.0]*len(TestTimes)
Correlation=[0.0]*len(TestTimes)
for x in range(0,len(TestTimes)):
    TimeStart=int(TestTimes[x])
    Tau,Corre,err=TDECorrelate(arrayA,arrayB,TimeStart,TimeBlock,StepSize,StepNumber)
    dt[x]=Tau
    Correlation[x]=Corre
    errors[x]=err
    print "T{} {}+/-{} {}%Correlation".format(x+1,Tau,err,Corre*100)
"""
