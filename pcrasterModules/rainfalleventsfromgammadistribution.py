from generalfunctions import *
from pcraster import *
from pcraster.framework import *
import sys
import random

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

class RainfallEventsFromGammaDistribution:
  def __init__(self, probabilityOfARainstorm, durationOfRainstorm, expectedRainfallIntensity, gammaShapeParameter, \
               timeStepsToReport,setOfVariablesToReport):
    '''Generates a random rainstorm event
    probabilityOfARainstorm - probability of occurence of a rainstorm when getRainstorm is called
    durationOfRainstorm - hours
    expectedRainfallIntensity - metre per hour
    '''
    self.probabilityOfARainstorm=probabilityOfARainstorm
    self.durationOfRainstorm=durationOfRainstorm
    self.expectedRainfallIntensity=expectedRainfallIntensity
    self.gammaShapeParameter=gammaShapeParameter

    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport

  def report(self,sample,timestep):
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      probRain=scalar(self.probabilityOfARainstorm)
      self.variablesToReport = {
                                'Pf': self.rainfallFlux,
                                'Pp': probRain
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = { }

    if timestep in self.timeStepsToReport:
      for variable in self.variablesToReport:
        report(self.variablesToReport[variable],generateNameST(variable, sample, timestep))

  def setProbabilityOfARainstorm(self,probabilityOfARainstorm):
    self.probabilityOfARainstorm=probabilityOfARainstorm

  def getRainstorm(self):
    import generalfunctions
    '''Create a realization of a rainstorm
    isRaining - Python integer true or false
    '''
    self.isRaining=random.uniform(0,1) < self.probabilityOfARainstorm 
    if self.isRaining:
      self.rainfallFlux=scalar(self.expectedRainfallIntensity*generalfunctions.mapgamma(self.gammaShapeParameter))
      self.rainfallAmount=scalar(self.rainfallFlux*self.durationOfRainstorm)
    else:
      self.rainfallFlux=scalar(0)
      self.rainfallAmount=scalar(0)  
    return self.isRaining, self.rainfallFlux, self.rainfallAmount

