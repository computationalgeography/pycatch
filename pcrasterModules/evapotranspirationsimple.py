from pcraster import *
from pcraster.framework import *

# Notes:
# time step duration in h
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a waterslice over a whole cell)
# if unit cannot be derived in this way (e.g. flux/fluxAmount/store), unit is indicated
# inputs of function is PCRaster type, inside function Python types are used

#setclone('clone.map')

class EvapotranspirationSimple:
  def __init__(self, timeStepDuration, beta, maximumEvapotranspirationFlux, timeStepsToReport, setOfVariablesToReport):
    '''
    beta, half saturation constant relating transpiration to vegetation density (kg/m2), typical 7 kg/m2
          where biomass is typically between zero and about 12
    maximumEvapotranspirationFlux, maximum evapotranspiration flux m/h
          500-1000 mm/year (niko et al, rhemedus site)
    '''

    self.timeStepDuration=timeStepDuration
    self.beta=beta
    self.maximumEvapotranspirationFlux=maximumEvapotranspirationFlux
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport

  def report(self,sample,timestep):
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                            'Ep': self.potentialEvapotranspirationFlux
                            #'Ea': self.actualEvapotranspirationFlux
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = {}

    if timestep in self.timeStepsToReport:
      for variable in self.variablesToReport:
        report(self.variablesToReport[variable],generateNameST(variable, sample, timestep))

  def potentialEvapotranspiration(self,fWaterPotential,biomass):
    self.potentialEvapotranspirationFlux=fWaterPotential*self.maximumEvapotranspirationFlux* \
                                    (biomass/(biomass+self.beta))
    return self.potentialEvapotranspirationFlux

  def actualEvapotranspiration(self,actualEvapotranspirationFlux):
    '''
    empty function to get actual evapo inside this class
    '''
    self.actualEvapotranspirationFlux=actualEvapotranspirationFlux
