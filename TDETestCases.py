"""
Compilation of functions used to make test cases
"""
import numpy as np
import random
pi=np.pi
Debugger=0
def PlaneWavePacket(Amp,k,omega,theta,sigma,x,y,SNR,length):
    Vx=(omega/k)*np.cos(theta)
    Vy=(omega/k)*np.sin(theta)
    kx=k*np.cos(theta)
    ky=k*np.sin(theta)
    t=np.arange(length)-int(length/2)
    sigmaPart=2*np.square(sigma)
    sine=Amp*np.cos((kx*x)+(ky*y)-(omega*t))
    #sine=(len(t)*[1.0]) (Uncomment to get wavepacket. Comment to just get gaussian)
    packet=np.exp((-np.square(x-(Vx*t))/sigmaPart))*np.exp((-np.square(y-(Vy*t))/sigmaPart))
    wavePacket=sine*packet
    maxAmp=max(wavePacket)
    noise=NoiseMaker(length,maxAmp,SNR)
    wavePacket=wavePacket+noise
    """
    t=np.arange(length)-int(length/2)
    tshiftX=t-shiftX
    tshiftY=t-shiftY
    sine=Amp*np.sin(2*pi*fX*tshiftX)*np.sin(2*pi*fY*tshiftY)
    Norm=1
    packet=Norm*np.exp(-np.square(tshiftX)/(2*np.square(sigmaX)))*np.exp(-np.square(tshiftY)/(2*np.square(sigmaY)))
    wavePacket=sine*packet
    maxAmp=max(wavePacket)
    noise=NoiseMaker(length,maxAmp,SNR)
    if (Debugger!=0):
        print (len(noise)-len(wavePacket))
    wavePacket=wavePacket+noise
    """
    return wavePacket

def ThreePointGenerator():
    Amp=100.0
    Base=20.0
    V=Base/50.0
    x=Base
    #print "point 1:({},0)".format(x)
    y=Base
    #print "point 2:(0,{})".format(x)
    theta=pi/4
    f=100.0
    omega=f
    sigma=10.0
    k=omega/V
    Vx=V*np.cos(theta)
    Vy=V*np.sin(theta)
    
    #print "Vx={}".format(Vx)
    #print "Vy={}".format(Vy)
    SNR=10.0
    length=1000
    t=np.arange(length)-(int(length/2))
    originPacket=PlaneWavePacket(Amp,k,omega,theta,sigma,0,0,SNR,length)
    dxPointPacket=PlaneWavePacket(Amp,k,omega,theta,sigma,x,0,SNR,length)
    dyPointPacket=PlaneWavePacket(Amp,k,omega,theta,sigma,0,y,SNR,length)
    return originPacket,dxPointPacket,dyPointPacket,Vx,Vy,x,y,t

def NoiseMaker(length,Amp,SNR):
    noise=np.array([0.0]*length)
    for x in range(0,len(noise)):
        noise[x]=random.gauss(0,Amp/SNR)
    return noise

def ThreePoint_NRunGenerator(N):
    #Creates an Nxlength matrix for 3 data points
    #Meant to simulate a concatination of N runs for 3 points
    for i in range(0,N):
        origin,dxPoint,dyPoint,Vx,Vy,x,y,t=ThreePointGenerator()
        if (i==0):
            originMatrix=origin
            xMatrix=dxPoint
            yMatrix=dyPoint
        else:
            originMatrix=np.vstack((originMatrix,origin))
            xMatrix=np.vstack((xMatrix,dxPoint))
            yMatrix=np.vstack((yMatrix,dyPoint))
    return originMatrix,xMatrix,yMatrix,Vx,Vy,x,y,t

class pixel:
    def __init__(self,xCoor,yCoor,timeData):
        #Initial Conditions
        self.xCoor=xCoor
        self.yCoor=yCoor
        self.timeData=timeData
        self.averageData=np.mean(timeData,axis=0)
        if (isinstance(self.averageData,(list,tuple,np.ndarray))==0):
            self.averageData=timeData
        
        #dt and Correlation
        self.dt=0.0
        self.errorDT=0.0
        self.Correlation=0.0
        self.errorCorrelation=0.0
        
        #Velocity
        self.Vx=0.0
        self.Vy=0.0
        self.errorVX=0.0
        self.errorVY=0.0

    def dtAndCorrelation(self,dt,errorDT,Correlation,errorCorrelation):
        self.dt=dt
        self.errorDT=errorDT
        self.Correlation=Correlation
        self.errorCorrelation=errorCorrelation

    def velocityRecorder(self,Vx,Vy,errorVX,errorVY):
        self.Vx=Vx
        self.Vy=Vy
        self.errorVX=errorVX
        self.errorVY=errorVY

    def Printer(self):
        x=self.xCoor
        y=self.yCoor
        dt=self.dt
        errorDT=self.errorDT
        Corre=self.Correlation
        eCorre=self.errorCorrelation
        Vx=self.Vx
        errorVX=self.errorVX
        Vy=self.Vy
        errorVY=self.errorVY
        print "Measurement at ({},{})".format(x,y)
        print "dt={}+/-{}".format(dt,errorDT)
        print "correlation={}+/-{}".format(Corre,eCorre)
        print "Vx={}+/-{}".format(Vx,errorVX)
        print "Vy={}+/-{}".format(Vy,errorVY)

    def dtAnalyzer(self,theoryDT,maxDT,maxErrorDT):
        theoryDiff=self.dt-theoryDT
        maxDiff=self.dt-maxDT
        if (theoryDT==0.0 and theoryDiff!=0):
            theoryDT=0.001
        if (maxDT==0.0 and maxDiff!=0):
            maxDT=0.001

        if (self.errorDT>=self.dt):
            print "Imprecise dt Measurement"
        if (self.Correlation<51):
            print "Low Correlation"
        if (self.errorCorrelation>=self.Correlation):
            print "Imprecise Correlation"
        
        if (theoryDiff==0.0):
            print "dt measurement is exact to theory"
        else:
            print "dt measurement is {}% from theoryDT".format(theoryDiff*(100.0/theoryDT))
            if (theoryDT<(self.dt-self.errorDT) or theoryDT>(self.dt+self.errorDT)):
                print "theoryDT is outside of errorbars"
                if (self.errorDT!=0):
                    print "{} errorbars from theoryDT".format(abs(theoryDiff)/self.errorDT)
        
        if (maxDiff==0.0):
            print "dt measurement is exact to maxChecker"
        else:
            print "dt measurement is {}% from maxDT".format((self.dt-maxDT)*(100.0/maxDT))
            if ((maxDT+maxErrorDT)<(self.dt-self.errorDT) or (maxDT-maxErrorDT)>(self.dt+self.errorDT)):
                print "maxDT is outside of errorbars"
                if (self.errorDT!=0):
                    print "{} errorbars from <maxDT>".format(abs(maxDiff)/self.errorDT)
    
    def velocityAnalyzer(self,theoryVX,maxVX,maxErrorVX,theoryVY,maxVY,maxErrorVY):
        theoryXDiff=self.Vx-theoryVX
        maxXDiff=self.Vx-maxVX
        theoryYDiff=self.Vy-theoryVY
        maxYDiff=self.Vy-maxVY
        if (theoryVX==0.0 and theoryXDiff!=0):
            theoryVX=0.001
        if (maxVX==0.0 and maxXDiff!=0):
            maxVX=0.001
        if (theoryVY==0.0 and theoryYDiff!=0):
            theoryVY=0.001
        if (maxVY==0.0 and maxYDiff!=0):
            maxVY=0.001
        if (self.errorVX>=self.Vx):
            print "Imprecise Vx Measurement"
        
        if (theoryXDiff==0.0):
            print "Vx measurement is exact to theory"
        else:
            print "Vx measurement is {}% from theoryVX".format((theoryXDiff)*(100.0/theoryVX))
            if (theoryVX<(self.Vx-self.errorVX) or theoryVX>(self.Vx+self.errorVX)):
                print "theoryVX is outside of errorbars"
                if (self.errorVX!=0):
                    print "{} errorbars from theoryVX".format(abs(theoryXDiff)/self.errorVX)
        
        if (maxXDiff==0.0):
            print "Vx measurement is exact to maxChecker"
        else:
            print "Vx measurement is {}% from maxVX".format((maxXDiff)*(100.0/maxVX))
            if ((maxVX+maxErrorVX)<(self.Vx-self.errorVX) or (maxVX-maxErrorVX)>(self.Vx+self.errorVX)):
                print "maxVX is outside of errorbars"
                if (self.errorVX!=0):
                    print "{} errorbars from <maxVX>".format(abs(maxXDiff)/self.errorVX)
        
        if (self.errorVY>=self.Vy):
            print "Imprecise Vy Measurement"

        if (theoryYDiff==0.0):
            print "Vy measurement is exact to theory"
        else:
            print "Vy measurement is {}% from theoryVY".format((theoryYDiff)*(100.0/theoryVY))
            if (theoryVY<(self.Vy-self.errorVY) or theoryVY>(self.Vy+self.errorVY)):
                print "theoryVY is outside of errorbars"
                if (self.errorVY!=0):
                    print "{} errorbars from theoryVY".format(abs(theoryYDiff)/self.errorVY)
        
        if (maxYDiff==0.0):
            print "Vy measurement is exact to maxChecker"
        else:
            print "Vy measurement is {}% from maxVY".format((maxYDiff)*(100.0/maxVY))
            if ((maxVY+maxErrorVY)<(self.Vy-self.errorVY) or (maxVY-maxErrorVY)>(self.Vy+self.errorVY)):
                print "maxVY is outside of errorbars"
                if (self.errorVY!=0):
                    print "{} errorbars from <maxVY>".format(abs(maxYDiff)/self.errorVY)
