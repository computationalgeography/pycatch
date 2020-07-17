#from generalfunctions import *
#import sys
#import datetime
import math
from PCRaster.Framework import *
from PCRaster import *

# notes
# time step duration in h
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a waterslice over a whole cell)
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

# note that south should be bottom of map!!
# note that units of x coor should be equal to units of elevation on digital elevation map

class PrecipitationFixedRainstorm:
  def __init__(self, endTimeStepOfRainstorm, precipitationFluxDuringRainstorm):
    self.endTimeStepOfRainstorm=endTimeStepOfRainstorm
    self.precipitationFluxDuringRainstorm=precipitationFluxDuringRainstorm


  def getPrecipitationFlux(self,currentTimeStep):
    if currentTimeStep < self.endTimeStepOfRainstorm:
      self.precipitationFlux=scalar(self.precipitationFluxDuringRainstorm)
    else:
      self.precipitationFlux=scalar(0)
    return self.precipitationFlux

  def report(self,sample,timestep):
    report(self.precipitationFlux,generateNameST('Po', sample, timestep))

  def printit(self):
    print 'solar altitude: ', self.solarAltitudeDegrees, 'solar azimuth: ', self.solarAzimuthDegrees
    print 'solar azimuth radians, South is pi: ', self.solarAzimuthRadiansConverted
