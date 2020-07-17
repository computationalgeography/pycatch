from pcraster import * 
import sys
from pcraster.framework import *
import component

# notes
# time step duration in h
# heads in m
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a waterslice over a whole cell)
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

class InfiltrationGreenAndAmpt(component.Component):
  def __init__(self, porosityFraction, initialMoistureContentFraction, saturatedConductivityFlux, 
               suctionHead, timeStepDuration,timeStepsToReport,setOfVariablesToReport):

    # init only for suspend and resume in filter
    self.potentialInfiltrationFlux=scalar(0)
    self.potentialInfiltrationFluxMetrePerTimeStep=scalar(0)
    self.variablesToReport={}

    # real inits
    self.porosityFraction=scalar(porosityFraction)
    self.initialMoistureContentFraction=scalar(initialMoistureContentFraction)
    self.timeStepDuration=scalar(timeStepDuration)
    self.saturatedConductivityFlux=saturatedConductivityFlux
    self.saturatedConductivityMetrePerTimeStep=scalar(self.saturatedConductivityFlux)*self.timeStepDuration
    self.suctionHead=scalar(suctionHead) # negative
    self.availablePoreSpace=self.porosityFraction-self.initialMoistureContentFraction
    self.verySmallValue=scalar(0.0000000000000000000001)
    self.actualInfiltrationFlux=self.verySmallValue
    self.store=self.verySmallValue  # i.e., cumulative infiltration

    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport

    # budget check
    self.initialStore=scalar(0.0)
    self.actualAdditionCum=scalar(0.0)

  def reportAsMaps(self,sample,timestep):
    # reports
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                               'Ii': self.actualInfiltrationFlux,
                               'Ij': self.potentialInfiltrationFlux,
                               'Is': self.store,
                               'Iks': self.saturatedConductivityFlux
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = { }
    self.reportMaps(sample,timestep)

  def updateVariablesAsNumpyToReport(self):
    self.variablesAsNumpyToReport = {
                                    }

  def reportAsNumpyOneFile(self,locations,sample,timestep,endTimeStep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpyOneFilePerRealization(locations,sample,timestep,endTimeStep)

  def reportAsNumpyMultipleFiles(self,locations,sample,timestep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpy(locations,sample,timestep)


  def amountToFlux(self,amount):
    flux=amount/self.timeStepDuration
    return flux

  def fluxToAmount(self,flux):
    amount=flux * self.timeStepDuration
    return amount

  def potentialInfiltrationFluxFunction(self):
    # potential infiltration per timestep ('rate', m/timestep)
    self.potentialInfiltrationFluxMetrePerTimeStep = \
         scalar( \
         self.saturatedConductivityMetrePerTimeStep * \
         ( ((scalar(0.0)-self.suctionHead)*self.availablePoreSpace+self.store) /  \
         self.store ))
    self.potentialInfiltrationFlux=self.potentialInfiltrationFluxMetrePerTimeStep/self.timeStepDuration
    return self.potentialInfiltrationFlux


  def update(self,availableForInfiltrationFlux):
    self.actualInfiltrationFlux=min(availableForInfiltrationFlux,self.potentialInfiltrationFlux)
    self.store=self.store+self.fluxToAmount(self.actualInfiltrationFlux)
    return self.actualInfiltrationFlux

  def printit(self, row, column):
    printCellValue(self, self.rainFlux, 'rainfall flux', 'm/h', row, column)

  def budgetCheck(self, sample, timestep):
    # NOTE this is only valid if addition,subtraction are invoked ONCE EACH TIME STEP
    # NOTE use of maptotal, in case of ldd not covering whole area, absolute values may not be 
    # comparable with other budgets.
    self.actualAdditionCum=scalar(self.actualAdditionCum+self.actualInfiltrationFlux*self.timeStepDuration)
    self.increaseInStore=self.store-self.initialStore
    budget=maptotal(self.actualAdditionCum-self.increaseInStore)
    report(budget,generateNameST('B-inf', sample, timestep))
    return self.increaseInStore
