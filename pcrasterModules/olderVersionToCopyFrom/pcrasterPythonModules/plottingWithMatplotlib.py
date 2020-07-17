import string
import os, shutil, numpy
import generalfunctions
import matplotlib.pyplot as plt
from PCRaster import *
from PCRaster.Framework import *
import aggregationfunctionsDK


sampleNumbers=range(1,25,1)
timeSteps=range(1,5,1)
row=1
col=1
data = aggregationfunctionsDK.mySelectTArrayForSamples("jan", sampleNumbers, timeSteps, row, col)

# plot realizations as timeseries
i = 0
while i < len(sampleNumbers):
  plt.plot(data[:,0],data[:,i+1])
  i = i + 1
plt.xlabel("time")
plt.ylabel("semivariance")
plt.gcf().clear()

# plot histogram of timestep
timestep=5
theTimestep = data[:,0] == timestep
# select the proper row
dataValues=numpy.compress(theTimestep,data,0)
#plt.hist(dataValues[0,:])
#plt.show()
print dataValues


arrayLower=aggregationfunctionsDK.mySelectTArrayNewNameFormatPercentilesWithTimesteps("jan", timeSteps, row, col,"0.1")
arrayHigher=aggregationfunctionsDK.mySelectTArrayNewNameFormatPercentilesWithTimesteps("jan", timeSteps, row, col,"0.9")
arrayMedian=aggregationfunctionsDK.mySelectTArrayNewNameFormatPercentilesWithTimesteps("jan", timeSteps, row, col,"0.5")

# plot median and conf interval
plt.fill_between(arrayLower[:,0],arrayLower[:,1], arrayHigher[:,1],facecolor="gray",linewidth=0)
plt.plot(arrayMedian[:,0],arrayMedian[:,1],linewidth=1,color="black")
plt.xlabel("time")
plt.ylabel("semivariance")
plt.gcf().clear()
#plt.show()

numpy.savetxt("test.txt",data)
#test=numpy.loadtxt("test.txt")
#print test
