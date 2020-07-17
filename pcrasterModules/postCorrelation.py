from pcraster import *
import aggregationfunctions

# please check aggregationfunctions.correlation and its included
# functions b4 using
nrSamples=200
nrTimeSteps=959
location=boolean('c02.map')  # location where correlation is calculated
dependentName='eh'  # variable (filename suffix)
independentName='th0_'   # variable (filename suffix)
locationName='test'
sampleNumbers=range(1,nrSamples+1,1)
#timeSteps=range(1,nrTimeSteps+1,1)
timeSteps=range(329,329+1,1)

aggregationfunctions.correlation(location, dependentName, independentName, locationName, sampleNumbers, timeSteps)
