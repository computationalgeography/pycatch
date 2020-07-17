# -*- coding: utf-8 -*-
from pcraster import * 
import sys, generalfunctions
from pcraster.framework import *
import component

# class to report random variables, mainly
# they are stored here and also resumed at a filter moment
# in addition they are stored at the last time step

class RandomParameters(component.Component):
  def __init__(self,timeStepsToReport, \
                    setOfVariablesToReport,  \
                    maximumInterceptionCapacityPerLAI, \
                    ksat,
                    regolithThicknessHomogeneous, \
                    saturatedConductivityMetrePerDay, \
                    multiplierMaxStomatalConductance):
    
    # init for variables required for filtering only (suspend, resume)
    # ..
    self.variablesAsNumpyToReport={}
    self.variablesToReport={}
 
    # real inits
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport
    self.maximumInterceptionCapacityPerLAI=maximumInterceptionCapacityPerLAI
    self.ksat=ksat
    self.regolithThicknessHomogeneous=regolithThicknessHomogeneous
    self.saturatedConductivityMetrePerDay=saturatedConductivityMetrePerDay
    self.multiplierMaxStomatalConductance=multiplierMaxStomatalConductance


  def reportAsMaps(self,sample,timestep):
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                                 'RPic': self.maximumInterceptionCapacityPerLAI,
                                 'RPks': self.ksat,
                                 'RPrt': self.regolithThicknessHomogeneous,
                                 'RPsc': self.saturatedConductivityMetrePerDay,
                                 'RPmm': self.multiplierMaxStomatalConductance
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = {
                                 'RPic': self.maximumInterceptionCapacityPerLAI,
                                 'RPks': self.ksat,
                                 'RPrt': self.regolithThicknessHomogeneous,
                                 'RPsc': self.saturatedConductivityMetrePerDay,
                                 'RPmm': self.multiplierMaxStomatalConductance
                                 }

    if timestep in self.timeStepsToReport:
      for variable in self.variablesToReport:
        report(self.variablesToReport[variable],generateNameST(variable, sample, timestep))

  def reportAsNumpy(self,locations,sample,timestep):
    self.variablesAsNumpyToReport = {
                                 'RPic': self.maximumInterceptionCapacityPerLAI,
                                 'RPks': self.ksat,
                                 'RPrt': self.regolithThicknessHomogeneous,
                                 'RPsc': self.saturatedConductivityMetrePerDay,
                                 'RPmm': self.multiplierMaxStomatalConductance
                                    }
    import generalfunctions
    for variable in self.variablesAsNumpyToReport:
      generalfunctions.reportLocationsAsNumpyArray(locations,self.variablesAsNumpyToReport[variable], \
                       variable,sample,timestep)

  def reportAsNumpyPostmcloop(self,samples,timeSteps):
    self.variablesAsNumpyToReport = {
                                 'RPic': self.maximumInterceptionCapacityPerLAI,
                                 'RPks': self.ksat,
                                 'RPrt': self.regolithThicknessHomogeneous,
                                 'RPsc': self.saturatedConductivityMetrePerDay,
                                 'RPmm': self.multiplierMaxStomatalConductance
                                    }
    import generalfunctions
    for variable in self.variablesAsNumpyToReport:
      aVariable = generalfunctions.openSamplesAndTimestepsAsNumpyArraysAsNumpyArray(variable,samples,timeSteps)
      numpy.save(variable,aVariable)

  def reportAtLastTimeStep(self,sample,timestep,nrofsteps):
    self.variablesToReport = {
                                 'RPic': self.maximumInterceptionCapacityPerLAI,
                                 'RPks': self.ksat,
                                 'RPrt': self.regolithThicknessHomogeneous,
                                 'RPsc': self.saturatedConductivityMetrePerDay,
                                 'RPmm': self.multiplierMaxStomatalConductance
                                }
    import generalfunctions
    for variable in self.variablesToReport:
      generalfunctions.reportAVariableAtTheLastTimeStep(timestep,sample,nrofsteps, self.variablesToReport[variable],variable)
