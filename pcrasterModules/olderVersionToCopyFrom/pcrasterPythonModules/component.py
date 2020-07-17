from PCRaster import *
from PCRaster.Framework import *

# notes
# time step duration in h
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a waterslice over a whole cell)
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

class Component:
  def __init__(self):
    pass

  def report(self, sample, timestep):
    for variable in self.variablesToReport:
      if timestep in self.timeStepsToReport:
        report(self.variablesToReport[variable],generateNameST(variable, sample, timestep))
