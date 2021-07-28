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

class InfiltrationOnlyKsat(component.Component):
  def __init__(self, saturatedConductivityFlux, bareSoilSaturatedConductivityFlux, \
               maxSaturatedConductivityFluxFromVegetation,biomassHalfSaturation, timeStepDuration,timeStepsToReport,setOfVariablesToReport):

    # init only for suspend and resume in filter
    self.variablesToReport={}

    # real inits
    self.timeStepDuration=scalar(timeStepDuration)
    self.saturatedConductivityFlux=saturatedConductivityFlux
    self.bareSoilSaturatedConductivityFlux=bareSoilSaturatedConductivityFlux
    self.maxSaturatedConductivityFluxFromVegetation=maxSaturatedConductivityFluxFromVegetation
    self.biomassHalfSaturation=biomassHalfSaturation
    self.verySmallValue=scalar(0.0000000000000000000001)
    self.actualInfiltrationFlux=self.verySmallValue
    self.store=self.verySmallValue  # i.e., cumulative infiltration

    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport

    # budget check
    self.initialStore=scalar(0.0)
    self.actualAdditionCum=scalar(0.0)

    self.output_mapping = {
                           'Ii': self.actualInfiltrationFlux,
                           'Is': self.store,
                           'Iks': self.saturatedConductivityFlux
                          }

  def reportAsMaps(self, sample, timestep):
    self.variablesToReport = self.rasters_to_report(self.setOfVariablesToReport)
    self.reportMaps(sample, timestep)

  def amountToFlux(self,amount):
    flux=amount/self.timeStepDuration
    return flux

  def fluxToAmount(self,flux):
    amount=flux * self.timeStepDuration
    return amount

  def potentialInfiltrationFluxFunction(self):
    return self.saturatedConductivityFlux

  def update(self,availableForInfiltrationFlux):
    self.actualInfiltrationFlux=min(availableForInfiltrationFlux,self.saturatedConductivityFlux)
    self.store=self.store+self.fluxToAmount(self.actualInfiltrationFlux)
    return self.actualInfiltrationFlux

  def setSaturatedConductivityFluxAsFunctionOfBiomass(self,biomass):
    saturatedConductivityFluxFromVegetation=self.maxSaturatedConductivityFluxFromVegetation*(biomass/(self.biomassHalfSaturation+biomass))
    self.saturatedConductivityFlux=saturatedConductivityFluxFromVegetation+self.bareSoilSaturatedConductivityFlux

  def budgetCheck(self, sample, timestep):
    # NOTE this is only valid if addition,subtraction are invoked ONCE EACH TIME STEP
    # NOTE use of maptotal, in case of ldd not covering whole area, absolute values may not be 
    # comparable with other budgets.
    self.actualAdditionCum=scalar(self.actualAdditionCum+self.actualInfiltrationFlux*self.timeStepDuration)
    self.increaseInStore=self.store-self.initialStore
    budget=maptotal(self.actualAdditionCum-self.increaseInStore)
    report(budget,generateNameST('B-inf', sample, timestep))
    return self.increaseInStore
