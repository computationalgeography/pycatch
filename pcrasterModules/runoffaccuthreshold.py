import pcraster as pcr
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
# the module has to update self.store by itself!
# also, the module has to update cumulative discharge by itself

class RunoffAccuthreshold(component.Component):
  def __init__(self, ldd, timeStepDuration, timeStepsToReport, setOfVariablesToReport):

    # init for supsend and resume in filtering
    self.RunoffCubicMetrePerHour = pcr.scalar(0)
    self.RunoffCubicMeterPerHourMovingAverage = pcr.scalar(0)
    self.actualAbstractionFlux = pcr.scalar(0)
    self.variablesToReport = {}
    self.variablesAsNumpyToReport = {}

    # real init
    self.ldd = ldd
    self.timeStepDuration = pcr.scalar(timeStepDuration)
    self.cellArea = pcr.cellarea()
    self.cumulativeDischargeCubicMetres = pcr.scalar(0.0)

    self.timeStepsToReport = timeStepsToReport
    self.setOfVariablesToReport = setOfVariablesToReport

  def reportAsMaps(self, sample, timestep):
    '''
    map report
    '''
    self.output_mapping = {
                            'Rq': self.RunoffCubicMetrePerHour,
                            'Rqs': self.RunoffCubicMeterPerHourMovingAverage
                          }
    self.variablesToReport = self.rasters_to_report(self.setOfVariablesToReport)
    self.reportMaps(sample, timestep)

  def updateVariablesAsNumpyToReport(self):
    self.variablesAsNumpyToReport = {
                                'Rq': self.RunoffCubicMetrePerHour,
                                'Rqs': self.RunoffCubicMeterPerHourMovingAverage
                                    }

  def reportAsNumpyOneFile(self, locations, sample, timestep, endTimeStep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpyOneFilePerRealization(locations, sample, timestep, endTimeStep)

  def reportAsNumpyMultipleFiles(self, locations, sample, timestep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpy(locations, sample, timestep)




  def setSurfaceProperties(self, ldd):
    self.ldd = ldd

  def perHourToPerSecond(self, aScalar):
    perSecond = aScalar / pcr.scalar(3600)
    return perSecond

  def perSecondToPerHour(self, aScalar):
    perHour = aScalar * pcr.scalar(3600)
    return perHour

  def perHourToPerTimeStep(self, aScalar):
    perTimeStep = aScalar * self.timeStepDuration
    return perTimeStep

  def perSecondToPerTimeStep(self, aScalar):
    perHour = self.perSecondToPerHour(aScalar)
    perTimeStep = self.perHourToPerTimeStep(perHour)
    return perTimeStep

  def updateCumulativeDischarge(self):
    self.cumulativeDischargeCubicMetres = self.cumulativeDischargeCubicMetres + \
                                        self.perHourToPerTimeStep(self.RunoffCubicMetrePerHour)

  def update(self, rainFlux, potentialAbstractionFlux):
    # moving average calculation
    runoffPreviousTimeStep = self.RunoffCubicMetrePerHour
    # runoff calculation
    self.RunoffCubicMetrePerHour = pcr.accuthresholdflux(self.ldd, rainFlux, potentialAbstractionFlux) * self.cellArea
    self.actualAbstractionFlux = pcr.accuthresholdstate(self.ldd, rainFlux, potentialAbstractionFlux)
    self.updateCumulativeDischarge()
    # moving average calculation
    self.RunoffCubicMeterPerHourMovingAverage = (runoffPreviousTimeStep + self.RunoffCubicMetrePerHour) / pcr.scalar(2.0)
    return self.actualAbstractionFlux, self.RunoffCubicMetrePerHour

  def budgetCheck(self):
    return self.cumulativeDischargeCubicMetres

# test
# setclone('clone.map')
# ldd='ldd.map'
# timestepDuration=2.0
#
#d_runoffAccuthreshold=RunoffAccuthreshold(ldd, timestepDuration)
#
# actualAbstractionFlux=d_runoffAccuthreshold.update(0.003,0.005)
#
# report(actualAbstractionFlux,'aaf.map')
# report(d_runoffAccuthreshold.RunoffCubicMetrePerHour,'q.map')
