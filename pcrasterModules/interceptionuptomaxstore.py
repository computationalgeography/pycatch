import pcraster as pcr
import pcraster.framework as pcrfw
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
  def __init__(self, ldd, initialStore, maximumStore, gapFraction, timeStepDuration, timeStepsToReport, setOfVariablesToReport):

    # init only to run supsend and resume in filtering
    self.variablesToReport = {}
    self.variablesAsNumpyToReport = {}

    # real inits
    self.ldd = ldd
    self.maximumStore = pcr.scalar(maximumStore)  # maximum storage over whole cell (not just 1-gapFraction)
    self.gapFraction = pcr.scalar(gapFraction)
    self.initialStore = pcr.scalar(initialStore)
    self.store = pcr.scalar(self.initialStore)
    self.timeStepDuration = pcr.scalar(timeStepDuration)
    self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour = pcr.scalar(0)

    # budget checks
    self.actualAdditionFlux = pcr.scalar(0.0)
    self.actualAbstractionFlux = pcr.scalar(0.0)
    self.actualAdditionCum = pcr.scalar(0.0)
    self.actualAbstractionCum = pcr.scalar(0.0)
    self.timeStepsToReport = timeStepsToReport
    self.setOfVariablesToReport = setOfVariablesToReport

    self.output_mapping = {
                                'Vs': self.store,
                                'Vo': self.actualAbstractionFlux,
                                'Vi': self.actualAdditionFlux,
                                'Vgf': self.gapFraction,
                                'Vms': self.maximumStore,
                                'Vot': self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour
                          }


  def reportAsMaps(self, sample, timestep):
    self.variablesToReport = self.rasters_to_report(self.setOfVariablesToReport)
    self.reportMaps(sample, timestep)

  def updateVariablesAsNumpyToReport(self):
    self.variablesAsNumpyToReport = {
                                'Vot': self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour
                                    }

  def reportAsNumpyOneFile(self, locations, sample, timestep, endTimeStep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpyOneFilePerRealization(locations, sample, timestep, endTimeStep)

  def reportAsNumpyMultipleFiles(self, locations, sample, timestep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpy(locations, sample, timestep)


  def fluxToAmount(self, flux):
    fluxAmount = flux * self.timeStepDuration
    return fluxAmount

  def amountToFlux(self, amount):
    flux = amount / self.timeStepDuration
    return flux

  def calculateMaximumAdditionAmount(self):
    return pcr.max(0, self.maximumStore - self.store)

  def calculateMaximumAbstractionAmount(self):
    return pcr.max(0, self.store - 0.0)

  def gapFractionLoss(self, potentialAdditionFlux):
    canopyFraction = 1.0 - self.gapFraction
    potentialAdditionToCanopyFlux = potentialAdditionFlux * canopyFraction
    return potentialAdditionToCanopyFlux

  def addWater(self, potentialAdditionFlux):
    potentialAdditionToCanopyFlux = self.gapFractionLoss(potentialAdditionFlux)
    potentialAdditionAmount = self.fluxToAmount(potentialAdditionToCanopyFlux)
    maximumAdditionAmount = self.calculateMaximumAdditionAmount()
    actualAdditionAmount = pcr.min(potentialAdditionAmount, maximumAdditionAmount)
    self.actualAdditionFlux = self.amountToFlux(actualAdditionAmount)

    self.store = pcr.max(pcr.min(self.store + actualAdditionAmount, self.maximumStore), 0)
    return self.actualAdditionFlux

  def abstractWater(self, potentialAbstractionFlux):
    potentialAbstractionAmount = self.fluxToAmount(potentialAbstractionFlux)
    maximumAbstractionAmount = self.calculateMaximumAbstractionAmount()
    actualAbstractionAmount = pcr.min(potentialAbstractionAmount, maximumAbstractionAmount)
    self.actualAbstractionFlux = self.amountToFlux(actualAbstractionAmount)

    # for reporting
    self.totalActualAbstractionInUpstreamArea()

    self.store = pcr.max(pcr.min(self.store - actualAbstractionAmount, self.maximumStore), 0)
    return self.actualAbstractionFlux

  def totalActualAbstractionInUpstreamArea(self):
    self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour = pcr.accuflux(self.ldd, self.actualAbstractionFlux) * pcr.cellarea()

  def setGapFraction(self, gapFraction):
    self.gapFraction = pcr.scalar(gapFraction)

  def setMaximumStore(self, maximumStore):
    self.maximumStore = pcr.scalar(maximumStore)
    self.store = pcr.min(self.store, self.maximumStore)

#  def printit(self, row, column):
#    printCellValue(self, self.throughfallFlux, 'throughfall flux', 'm/h', row, column)
#    printCellValue(self, self.store, 'interception store', 'm', row, column)

  def budgetCheck(self, sample, timestep):
    # this should include setMaximumStore as this may result in throwing away of water
    # NOTE this is only valid if addition,subtraction and lateral flow are invoked ONCE EACH TIME STEP
    # NOTE use of maptotal, in case of ldd not covering whole area, absolute values may not be
    # comparable with other budgets.
    self.actualAdditionCum = self.actualAdditionCum + self.actualAdditionFlux * self.timeStepDuration
    self.actualAbstractionCum = self.actualAbstractionCum + self.actualAbstractionFlux * self.timeStepDuration
    self.increaseInStore = self.store - self.initialStore
    # budget=catchmenttotal(self.actualAdditionCum-self.actualAbstractionCum-self.increaseInStore)
    budget = pcr.maptotal(self.actualAdditionCum - self.actualAbstractionCum - self.increaseInStore)
    pcr.report(budget, pcrfw.generateNameST('B-int', sample, timestep))
    return self.increaseInStore

# setclone("clone.map")
# d_interceptionuptomaxstore=InterceptionUpToMaxStore(0.1,0.9,0.5,1.0,'test','test')
# jan=d_interceptionuptomaxstore.addWater(20.0)
# report(jan,"jan")
# report(d_interceptionuptomaxstore.store,"store1")
# piet=d_interceptionuptomaxstore.abstractWater(0.01)
# report(piet,"piet")
# report(d_interceptionuptomaxstore.store,"store2")
