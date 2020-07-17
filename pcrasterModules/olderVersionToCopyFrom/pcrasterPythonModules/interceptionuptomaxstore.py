from PCRaster import *
from PCRaster.Framework import *
import component

# notes
# time step duration in h
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a waterslice over a whole cell)
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

class InterceptionUpToMaxStore(component.Component):
  def __init__(self, initialStore, maximumStore, gapFraction, timeStepDuration,timeStepsToReport,setOfVariablesToReport):

    # init only to run supsend and resume in filtering
    self.variablesToReport={}

    # real inits
    self.maximumStore=scalar(maximumStore)  # maximum storage over whole cell (not just 1-gapFraction)
    self.gapFraction=scalar(gapFraction)
    self.initialStore=scalar(initialStore)
    self.store=scalar(self.initialStore)
    self.timeStepDuration=scalar(timeStepDuration)
    # budget checks
    self.actualAdditionFlux=scalar(0.0)
    self.actualAbstractionFlux=scalar(0.0)
    self.actualAdditionCum=scalar(0.0)
    self.actualAbstractionCum=scalar(0.0)
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport


  def report(self,sample,timestep):
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                                'Vs': self.store,
                                'Vo': self.actualAbstractionFlux,
                                'Vi': self.actualAdditionFlux,
                                'Vgf': self.gapFraction,
                                'Vms': self.maximumStore
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = {}

    if timestep in self.timeStepsToReport:
      for variable in self.variablesToReport:
        report(self.variablesToReport[variable],generateNameST(variable, sample, timestep))

  def fluxToAmount(self,flux):
    fluxAmount=flux * self.timeStepDuration
    return fluxAmount

  def amountToFlux(self,amount):
    flux=amount/self.timeStepDuration
    return flux

  def calculateMaximumAdditionAmount(self):
    return max(0,self.maximumStore-self.store)

  def calculateMaximumAbstractionAmount(self):
    return max(0,self.store-0.0)

  def gapFractionLoss(self,potentialAdditionFlux):
    canopyFraction=1.0-self.gapFraction
    potentialAdditionToCanopyFlux=potentialAdditionFlux*canopyFraction
    return potentialAdditionToCanopyFlux

  def addWater(self,potentialAdditionFlux):
    potentialAdditionToCanopyFlux=self.gapFractionLoss(potentialAdditionFlux)
    potentialAdditionAmount=self.fluxToAmount(potentialAdditionToCanopyFlux)
    maximumAdditionAmount=self.calculateMaximumAdditionAmount()
    actualAdditionAmount=min(potentialAdditionAmount,maximumAdditionAmount)
    self.actualAdditionFlux=self.amountToFlux(actualAdditionAmount)

    self.store=max(min(self.store+actualAdditionAmount,self.maximumStore),0)
    return self.actualAdditionFlux

  def abstractWater(self,potentialAbstractionFlux):
    potentialAbstractionAmount=self.fluxToAmount(potentialAbstractionFlux)
    maximumAbstractionAmount=self.calculateMaximumAbstractionAmount()
    actualAbstractionAmount=min(potentialAbstractionAmount,maximumAbstractionAmount)
    self.actualAbstractionFlux=self.amountToFlux(actualAbstractionAmount)
    
    self.store=max(min(self.store-actualAbstractionAmount,self.maximumStore),0)
    return self.actualAbstractionFlux

  def setGapFraction(self,gapFraction):
    self.gapFraction=scalar(gapFraction)

  def setMaximumStore(self,maximumStore):
    self.maximumStore=scalar(maximumStore) 
    self.store=min(self.store,self.maximumStore)

#  def printit(self, row, column):
#    printCellValue(self, self.throughfallFlux, 'throughfall flux', 'm/h', row, column)
#    printCellValue(self, self.store, 'interception store', 'm', row, column)

  def budgetCheck(self, sample, timestep):
    # this should include setMaximumStore as this may result in throwing away of water
    # NOTE this is only valid if addition,subtraction and lateral flow are invoked ONCE EACH TIME STEP
    # NOTE use of maptotal, in case of ldd not covering whole area, absolute values may not be 
    # comparable with other budgets.
    self.actualAdditionCum=self.actualAdditionCum+self.actualAdditionFlux*self.timeStepDuration
    self.actualAbstractionCum=self.actualAbstractionCum+self.actualAbstractionFlux*self.timeStepDuration
    self.increaseInStore=self.store-self.initialStore
    #budget=catchmenttotal(self.actualAdditionCum-self.actualAbstractionCum-self.increaseInStore)
    budget=maptotal(self.actualAdditionCum-self.actualAbstractionCum-self.increaseInStore)
    report(budget,generateNameST('B-int', sample, timestep))
    return self.increaseInStore

#setclone("clone.map")
#d_interceptionuptomaxstore=InterceptionUpToMaxStore(0.1,0.9,0.5,1.0,'test','test')
#jan=d_interceptionuptomaxstore.addWater(20.0)
#report(jan,"jan")
#report(d_interceptionuptomaxstore.store,"store1")
#piet=d_interceptionuptomaxstore.abstractWater(0.01)
#report(piet,"piet")
#report(d_interceptionuptomaxstore.store,"store2")
