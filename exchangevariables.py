# -*- coding: utf-8 -*-
from pcraster import * 
import sys
from pcraster.framework import *

class ExchangeVariables:
  def __init__(self,timeStepsToReport,setOfVariablesToReport):
    
    # init for variables required for filtering only (suspend, resume)
    self.cumulativePrecipitation=scalar(0)
    self.evapFromSoilMultiplier=scalar(0)
    self.regolithThicknessAParameter=scalar(0)
    self.regolithThicknessCParameter=scalar(0)
    self.upwardSeepageFlux=scalar(0)
    self.variablesToReport={}
 
    # real inits
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport

  def report(self,sample,timestep):
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                               'Xrc': self.regolithThicknessCParameter,
                               #'Xra': self.regolithThicknessAParameter,
                               #'Xmu': self.evapFromSoilMultiplier
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = {
                               'Xrc': self.regolithThicknessCParameter,
                               'Xra': self.regolithThicknessAParameter,
                                 }

    if timestep in self.timeStepsToReport:
      for variable in self.variablesToReport:
        report(self.variablesToReport[variable],generateNameST(variable, sample, timestep))
