from PCRaster import * 
import sys
from PCRaster.Framework import *

# notes
# time step duration in years
# levels in m
# vertical change in m/year, rise is positive
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

class Baselevel:
  def __init__(self, areaWhereBaselevelIsSet, initialLevel, baseLevelRise, timeStepDuration, \
               timeStepsToReport, setOfVariablesToReport):

    # real inits
    self.areaWhereBaselevelIsSet=areaWhereBaselevelIsSet
    self.initialLevel=initialLevel
    self.baseLevelRise=baseLevelRise
    self.timeStepDuration=timeStepDuration
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport
    self.baselevel=ifthen(self.areaWhereBaselevelIsSet,self.initialLevel)

  def report(self,sample,timestep):
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                               'Ll': self.baselevel
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = {
                               'Ll': self.baselevel
                                 }

    if timestep in self.timeStepsToReport:
      for variable in self.variablesToReport:
        report(self.variablesToReport[variable],generateNameST(variable, sample, timestep))

  def getBaselevel(self,timeStep):
    '''
    calculates baselevel, where it is not set, missing values
    '''
    riseSinceStart=self.baseLevelRise*timeStep*self.timeStepDuration
    self.baselevel=ifthen(self.areaWhereBaselevelIsSet,self.initialLevel+riseSinceStart)
    return self.baselevel


#  def budgetCheck(self, sample, timestep):
   # needs to take into account actualDepositionFlux
