from PCRaster import *
from PCRaster.Framework import *
import sys

# notes
# time step duration in h
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a waterslice over a whole cell)
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

class InterceptionFraction:
  def __init__(self, fractionIntercepted, timeStep):
    self.store=scalar(0)
    self.fractionIntercepted=scalar(fractionIntercepted)
    self.timeStep=timeStep

  def calculateThroughfallFlux(self, rainFlux):
    self.throughfallFlux=rainFlux*(1-self.fractionIntercepted)

  def calculateInterceptionStore(self, rainFlux):
    self.store = self.store + (rainFlux-self.throughfallFlux)*self.timeStep

  def update(self, rainFlux):
    self.calculateThroughfallFlux(rainFlux) 
    self.calculateInterceptionStore(rainFlux)
    return self.throughfallFlux

  def printit(self, row, column):
    printCellValue(self, self.throughfallFlux, 'throughfall flux', 'm/h', row, column)
    printCellValue(self, self.store, 'interception store', 'm', row, column)

  def report(self, sample, timestep):
    report(self.throughfallFlux,generateNameST('tf', sample, timestep))
    report(self.store,generateNameST('ints', sample, timestep))
