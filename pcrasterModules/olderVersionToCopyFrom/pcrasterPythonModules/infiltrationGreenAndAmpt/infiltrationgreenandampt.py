from PCRaster import *
import sys
from PCRaster.Framework import *
from generalfunctions import *

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

class InfiltrationGreenAndAmpt:
  def __init__(self, porosityFraction, initialMoistureContentFraction, saturatedConductivityFlux, 
               suctionHead, timeStepDuration):
    self.porosityFraction=scalar(porosityFraction)
    self.initialMoistureContentFraction=scalar(initialMoistureContentFraction)
    self.timeStepDuration=scalar(timeStepDuration)
    self.saturatedConductivityMetrePerTimeStep=scalar(saturatedConductivityFlux)*self.timeStepDuration
    self.suctionHead=suctionHead # negative
    self.availablePoreSpace=self.porosityFraction-self.initialMoistureContentFraction
    self.verySmallValue=scalar(0.0000000000000000000001)
    self.actualInfiltrationFlux=self.verySmallValue
    self.store=self.verySmallValue  # i.e., cumulative infiltration

    # budget check
    self.initialStore=0.0
    self.actualAdditionCum=0.0

  def amountToFlux(self,amount):
    flux=amount/self.timeStepDuration
    return flux

  def fluxToAmount(self,flux):
    amount=flux * self.timeStepDuration
    return amount

  def potentialInfiltrationFluxFunction(self):
    # potential infiltration per timestep ('rate', m/timestep)
    self.potentialInfiltrationFluxMetrePerTimeStep = \
         self.saturatedConductivityMetrePerTimeStep * \
         ( ((scalar(0.0)-self.suctionHead)*self.availablePoreSpace+self.store) /  \
         self.store )
    self.potentialInfiltrationFlux=self.potentialInfiltrationFluxMetrePerTimeStep/self.timeStepDuration
    return self.potentialInfiltrationFlux


  def update(self,availableForInfiltrationFlux):
    self.actualInfiltrationFlux=min(availableForInfiltrationFlux,self.potentialInfiltrationFlux)
    self.store=self.store+self.fluxToAmount(self.actualInfiltrationFlux)
    return self.actualInfiltrationFlux

  def report(self, sample, timestep):
    report(self.actualInfiltrationFlux,generateNameST('Ii', sample, timestep))
    report(mapaverage(self.actualInfiltrationFlux),generateNameST('Iia', sample, timestep))
    report(self.potentialInfiltrationFlux,generateNameST('Ij', sample, timestep))
    report(mapaverage(self.potentialInfiltrationFlux),generateNameST('Ija', sample, timestep))
    report(self.store,generateNameST('Is', sample, timestep))

  def printit(self, row, column):
    printCellValue(self, self.rainFlux, 'rainfall flux', 'm/h', row, column)

  def budgetCheck(self, sample, timestep):
    # NOTE this is only valid if addition,subtraction are invoked ONCE EACH TIME STEP
    # NOTE use of maptotal, in case of ldd not covering whole area, absolute values may not be 
    # comparable with other budgets.
    self.actualAdditionCum=self.actualAdditionCum+self.actualInfiltrationFlux*self.timeStepDuration
    self.increaseInStore=self.store-self.initialStore
    budget=maptotal(self.actualAdditionCum-self.increaseInStore)
    report(budget,generateNameST('B-inf', sample, timestep))
    return self.increaseInStore
