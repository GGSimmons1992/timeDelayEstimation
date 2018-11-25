#Error analysis calculation for 1/20/15 error analysis experiment
#More detail found in pages 28-30 of Gary Simmons' lab notebook
import numpy as np

errorA=np.array([428.8,96.0,51.51,75.33,287.0,34.9,9.35,25.9,34.5])
measurementA=np.array([380,130,160,220,120,50,10,20,30])
errorPercentA=np.mean(errorA/measurementA)
scoreA=errorPercentA/len(errorA)
print errorPercentA
print scoreA

errorB=np.array([40.4,63.1,91.5,147.7,116,202.6,180.6])
measurementB=np.array([20,60,10,210,70,220,120])
errorPercentB=np.mean(errorB/measurementB)
scoreB=errorPercentB/len(errorB)
print errorPercentB
print scoreB

errorC=np.array([261.9,126,63.7,285.1,59.8,175,168,93.3,81.8])
measurementC=np.array([160,70,10,360,50,30,100,50,10])
errorPercentC=np.mean(errorC/measurementC)
scoreC=errorPercentC/len(errorC)
print errorPercentC
print scoreC
