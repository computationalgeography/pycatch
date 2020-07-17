from pcraster import *
import sys, generalfunctions
from pcraster.framework import *

# notes
# time step duration in h
# vertical fluxes, variable name 'flux'
#   water: m/h (except where indicated)
#   vertical fluxes over a time step, variable name 'fluxAmount'
# amounts in storages, variable name 'store'
# energy in J/m2, 'Energy'
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

class BedrockWeathering:
  def __init__(self,weatheringRateBareBedrock,weatheringExponentParameter, timeStepsToReport,setOfVariablesToReport):
    '''
    weatheringRateBareBedrock in m/year
    weatheringExponentParameter in m-1
    '''
    self.weatheringRateBareBedrock=weatheringRateBareBedrock
    self.weatheringExponentParameter=weatheringExponentParameter
    self.weatheringMetrePerYear=spatial(scalar(0))
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport

  def report(self,sample,timestep):
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                                'Cwe': self.weatheringMetrePerYear,
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = {
                                'Cwe': self.weatheringMetrePerYear,
                                }

    if timestep in self.timeStepsToReport:
      for variable in self.variablesToReport:
        report(self.variablesToReport[variable],generateNameST(variable, sample, timestep))


  def weatheringRate(self,soilDepthMetres):
    self.weatheringMetrePerYear=self.weatheringRateBareBedrock * exp(-self.weatheringExponentParameter * soilDepthMetres)
    return self.weatheringMetrePerYear

  def steadyStateSoilDepth(self,baselevelFall):
    steadyStateSoilDepth=ln(self.weatheringRateBareBedrock/baselevelFall)/(self.weatheringExponentParameter)
    return steadyStateSoilDepth

### test
#dem='mdtpaz4.map'
#soilThick=ifthen(defined(dem),uniform(1)*10)
#d_bedrockweathering=BedrockWeathering(0.001,1)
#weathering=d_bedrockweathering.weatheringRate(soilThick)
#report(weathering,'weather.map')
#report(soilThick,'soilthick.map')

