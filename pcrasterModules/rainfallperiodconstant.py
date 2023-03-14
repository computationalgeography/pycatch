from pcraster import *
from generalfunctions import *
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
# class needs to update store by itself
# note that, by exception, rainfall store can be negative

class RainfallPeriodConstant:
  def __init__(self, constantRainFlux, startRainStorm, endRainStorm, timeStepDuration):
    self.constantRainFlux=scalar(constantRainFlux)
    self.startRainStorm=startRainStorm
    self.endRainStorm=endRainStorm
    self.timeStepDuration=timeStepDuration
    self.store=scalar(0)

  def perHourToPerTimeStep(self,aScalar):
    perTimeStep=aScalar*self.timeStepDuration
    return perTimeStep

  def updateStore(self):
    self.store=self.store-self.perHourToPerTimeStep(self.rainFlux)

  def update(self, currentTimeStep):
    raining=onePeriod(self,self.startRainStorm , self.endRainStorm, self.timeStepDuration, currentTimeStep)
    if raining:
      self.rainFlux=scalar(self.constantRainFlux)
    else:
      self.rainFlux=scalar(0)
    self.updateStore()
    return self.rainFlux

  def printit(self, row, column):
    printCellValue(self, self.rainFlux, 'rainfall flux', 'm/h', row, column)

  def report(self, sample, timestep):
    report(self.rainFlux,generateNameST('rain', sample, timestep))
