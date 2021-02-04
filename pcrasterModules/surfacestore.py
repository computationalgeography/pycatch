import pcraster as pcr
import pcraster.framework as pcrfw
import generalfunctions


# notes
# time step duration in h
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a waterslice over a whole cell)
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

maxErrorInFluxesForErrorMessages = 0.00001

class SurfaceStore:
  def __init__(self, initialStore, maxStore, timeStepDuration, timeStepsToReport, setOfVariablesToReport):

    # init only for supsend and resume in filtering
    self.changeAmount = pcr.scalar(0)
    self.variablesToReport = {}

    # real inits
    self.initialStore = initialStore
    self.store = pcr.scalar(self.initialStore)
    self.maxStore = maxStore
    self.changeFlux = pcr.scalar(0)
    self.timeStepDuration = pcr.scalar(timeStepDuration)
    self.timeStepsToReport = timeStepsToReport
    self.setOfVariablesToReport = setOfVariablesToReport

    # budget
    self.actualAdditionCum = pcr.scalar(0.0)


  def reportAsMaps(self, sample, timestep):
    # reports
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                                'Ss': self.store,
                                'Sc': self.changeFlux 
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = { }
    self.reportMaps(sample, timestep)

  def updateVariablesAsNumpyToReport(self):
    self.variablesAsNumpyToReport = {
                                    }

  def reportAsNumpyOneFile(self, locations, sample, timestep, endTimeStep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpyOneFilePerRealization(locations, sample, timestep, endTimeStep)

  def reportAsNumpyMultipleFiles(self, locations, sample, timestep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpy(locations, sample, timestep)



  def emptyIt(self):
    self.store = pcr.scalar(0.0)

  def getStore(self):
    return self.store

  def amountToFlux(self, amount):
    flux = amount / self.timeStepDuration
    return flux

  def fluxToAmount(self, flux):
    amount = flux * self.timeStepDuration
    return amount

  def potentialOutFlux(self):
    potentialOutFlux = self.amountToFlux(self.store)
    return potentialOutFlux

  def potentialToAmount(self):
    potentialToAmount = pcr.max(self.maxStore - self.store, pcr.scalar(0))
    return potentialToAmount

  def potentialToFlux(self):
    potentialToAmount = self.potentialToAmount()
    potentialToFlux = self.amountToFlux(potentialToAmount)
    return potentialToFlux

  def update(self, changeFlux):
    potentialToFlux = self.potentialToFlux()
    changeFluxGreaterThanPotentialToFlux = changeFlux >= potentialToFlux + maxErrorInFluxesForErrorMessages
    generalfunctions.printErrorMessageIfACellContainsTrue(changeFluxGreaterThanPotentialToFlux,
                                         'added more water to surface store than possible')
    self.changeAmount = self.fluxToAmount(changeFlux)
    self.changeFlux = self.amountToFlux(self.changeAmount)
    self.store = self.store + self.changeAmount

  def budgetCheck(self, sample, timestep):
    # NOTE this is only valid if addition,subtraction and lateral flow are invoked ONCE EACH TIME STEP
    # NOTE use of maptotal, in case of ldd not covering whole area, absolute values may not be 
    # comparable with other budgets.
    # note self.changeAmount can be positive or negative and thus self.actualAdditionCum, too
    self.actualAdditionCum = self.actualAdditionCum + self.changeAmount
    self.increaseInStore = self.store - self.initialStore
    # budget=catchmenttotal(self.actualAdditionCum-self.actualAbstractionCum-self.increaseInStore)
    budget = pcr.maptotal(self.actualAdditionCum - self.increaseInStore)
    pcr.report(budget, pcrfw.generateNameST('B-sur', sample, timestep))
    return self.increaseInStore

  # KDJ
  # def __repr__(self):
  #   return "%s %s %s ..." % (self.rainFlux, 'rainfall flux', 'm/h', row, column)

#  def printit(self, row, column):
#    printCellValue(self, self.rainFlux, 'rainfall flux', 'm/h', row, column)

# test
# setclone('clone.map')
# surfaceStore=ifthen('clone.map',pcr.scalar(0.2))
# maxStore=ifthen('clone.map',pcr.scalar(0.3))
#d_surfaceStore=SurfaceStore(surfaceStore, maxSurfaceStore, 1.0)
# d_surfaceStore.update(uniform('clone.map'))
