import pcraster as pcr
import pcraster.framework as pcrfw


class ExchangeVariables:
  def __init__(self, timeStepsToReport, setOfVariablesToReport):

    # init for variables required for filtering only (suspend, resume)
    self.cumulativePrecipitation = pcr.scalar(0)
    self.evapFromSoilMultiplier = pcr.scalar(0)
    self.regolithThicknessAParameter = pcr.scalar(0)
    self.regolithThicknessCParameter = pcr.scalar(0)
    self.upwardSeepageFlux = pcr.scalar(0)
    self.variablesToReport = {}

    # real inits
    self.timeStepsToReport = timeStepsToReport
    self.setOfVariablesToReport = setOfVariablesToReport

    self.output_mapping = {
                               'Xrc': self.regolithThicknessCParameter,
                               'Xra': self.regolithThicknessAParameter,
                               'Xmu': self.evapFromSoilMultiplier
                          }

  def report(self, sample, timestep):
    self.variablesToReport = self.rasters_to_report(self.setOfVariablesToReport)

    if timestep in self.timeStepsToReport:
      for variable in self.variablesToReport:
        pcr.report(self.variablesToReport[variable], pcrfw.generateNameST(variable, sample, timestep))
