import pcraster as pcr
import pcraster.framework as pcrfw
import sys, generalfunctions
import component

# class to report random variables, mainly
# they are stored here and also resumed at a filter moment
# in addition they are stored at the last time step

class RandomParameters(component.Component):
  def __init__(self, timeStepsToReport,
                    setOfVariablesToReport,
                    maximumInterceptionCapacityPerLAI,
                    ksat,
                    regolithThicknessHomogeneous,
                    saturatedConductivityMetrePerDay,
                    multiplierMaxStomatalConductance):

    # init for variables required for filtering only (suspend, resume)
    # ..
    self.variablesAsNumpyToReport = {}
    self.variablesToReport = {}

    # real inits
    self.timeStepsToReport = timeStepsToReport
    self.setOfVariablesToReport = setOfVariablesToReport
    self.maximumInterceptionCapacityPerLAI = maximumInterceptionCapacityPerLAI
    self.ksat = ksat
    self.regolithThicknessHomogeneous = regolithThicknessHomogeneous
    self.saturatedConductivityMetrePerDay = saturatedConductivityMetrePerDay
    self.multiplierMaxStomatalConductance = multiplierMaxStomatalConductance

    self.output_mapping = {
                                 'RPic': self.maximumInterceptionCapacityPerLAI,
                                 'RPks': self.ksat,
                                 'RPrt': self.regolithThicknessHomogeneous,
                                 'RPsc': self.saturatedConductivityMetrePerDay,
                                 'RPmm': self.multiplierMaxStomatalConductance
                          }


  def reportAsMaps(self, sample, timestep):
    self.variablesToReport = self.rasters_to_report(self.setOfVariablesToReport)

    if timestep in self.timeStepsToReport:
      for variable in self.variablesToReport:
        pcr.report(self.variablesToReport[variable], pcrfw.generateNameST(variable, sample, timestep))

  def reportAsNumpy(self, locations, sample, timestep):
    self.variablesAsNumpyToReport = {
                                 'RPic': self.maximumInterceptionCapacityPerLAI,
                                 'RPks': self.ksat,
                                 'RPrt': self.regolithThicknessHomogeneous,
                                 'RPsc': self.saturatedConductivityMetrePerDay,
                                 'RPmm': self.multiplierMaxStomatalConductance
                                    }
    import generalfunctions
    for variable in self.variablesAsNumpyToReport:
      generalfunctions.reportLocationsAsNumpyArray(locations, self.variablesAsNumpyToReport[variable],
                       variable, sample, timestep)

  def reportAsNumpyPostmcloop(self, samples, timeSteps):
    self.variablesAsNumpyToReport = {
                                 'RPic': self.maximumInterceptionCapacityPerLAI,
                                 'RPks': self.ksat,
                                 'RPrt': self.regolithThicknessHomogeneous,
                                 'RPsc': self.saturatedConductivityMetrePerDay,
                                 'RPmm': self.multiplierMaxStomatalConductance
                                    }
    import generalfunctions
    for variable in self.variablesAsNumpyToReport:
      aVariable = generalfunctions.openSamplesAndTimestepsAsNumpyArraysAsNumpyArray(variable, samples, timeSteps)
      numpy.save(variable, aVariable)

  def reportAtLastTimeStep(self, sample, timestep, nrofsteps):
    self.variablesToReport = {
                                 'RPic': self.maximumInterceptionCapacityPerLAI,
                                 'RPks': self.ksat,
                                 'RPrt': self.regolithThicknessHomogeneous,
                                 'RPsc': self.saturatedConductivityMetrePerDay,
                                 'RPmm': self.multiplierMaxStomatalConductance
                                }
    import generalfunctions
    for variable in self.variablesToReport:
      generalfunctions.reportAVariableAtTheLastTimeStep(timestep, sample, nrofsteps, self.variablesToReport[variable], variable)
