from PCRaster import *
from PCRaster.Framework import *

# notes
# time step duration in h
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a waterslice over a whole cell)
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

class InterceptionUpToMaxStore:
  def __init__(self, initialStore, maximumStore, timeStepDuration):
    self.maximumStore=scalar(maximumStore)
    self.initialStore=scalar(initialStore)
    self.store=scalar(self.initialStore)
    self.timeStepDuration=timeStepDuration
    # budget checks
    self.actualAdditionFlux=scalar(0.0)
    self.actualAbstractionFlux=scalar(0.0)
    self.actualAdditionCum=scalar(0.0)
    self.actualAbstractionCum=scalar(0.0)

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

  def addWater(self,potentialAdditionFlux):
    potentialAdditionAmount=self.fluxToAmount(potentialAdditionFlux)
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

#  def printit(self, row, column):
#    printCellValue(self, self.throughfallFlux, 'throughfall flux', 'm/h', row, column)
#    printCellValue(self, self.store, 'interception store', 'm', row, column)

  def report(self, sample, timestep):
    report(self.store,generateNameST('Vs', sample, timestep))
    report(self.actualAbstractionFlux,generateNameST('Vo', sample, timestep)) # note this is not throughfall
    report(self.actualAdditionFlux,generateNameST('Vi', sample, timestep))

  def budgetCheck(self, sample, timestep):
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
#d_interceptionuptomaxstore=InterceptionUpToMaxStore(0.1,0.9,1.0)
#jan=d_interceptionuptomaxstore.addWater(0.001)
#report(jan,"jan")
#report(d_interceptionuptomaxstore.store,"store1")
#piet=d_interceptionuptomaxstore.abstractWater(0.01)
#report(piet,"piet")
#report(d_interceptionuptomaxstore.store,"store2")
