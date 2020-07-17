#from generalfunctions import *
#import sys
#import datetime
import solar, math
from pcraster import *

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

class PrecipitationFixedRainstorms:
  def __init__(self, durationBetweenStartTimeOfRainstorms, durationOfARainstorm,precipitationFluxDuringRainstorm):
  """Generates precipitation for rainstorms with fixed rainfall intensity

  durationBetweenStartTimeOfRainstorms -- real duration between start of rainstorms, e.g. ten hours is timedelta(hours=10)
  durationOfARainstorm -- real duration of a rainstorm, e.g. two hours is timedelta(hours=2)
  precipitationFluxDuringRainstorm -- precipitation intensity during rainstorm (m/h)
  """


  # KDJ
  # def __repr__(self):
  #   return "%s %s %s ..." % (self.rainFlux, 'rainfall flux', 'm/h', row, column)

  def printit(self):
    print 'solar altitude: ', self.solarAltitudeDegrees, 'solar azimuth: ', self.solarAzimuthDegrees
    print 'solar azimuth radians, South is pi: ', self.solarAzimuthRadiansConverted
