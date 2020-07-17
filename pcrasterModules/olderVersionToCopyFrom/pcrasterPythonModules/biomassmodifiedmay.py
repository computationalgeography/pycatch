from PCRaster import *
import sys, generalfunctions
from PCRaster.Framework import *

# notes
# time step duration in hours
# vertical fluxes, variable name 'flux'
#   water: m/hour (except where indicated)
#   vertical fluxes over a time step, variable name 'fluxAmount'
# amounts in storages, variable name 'store'
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

class BiomassModifiedMay:
  def __init__(self, \
               biomass, \
               waterUseEfficiency, \
               maintenanceRate, \
               gamma, \
               alpha, \
               dispersion, \
               runoff, \
               sdOfNoise, \
               timeStepDuration,timeStepsToReport,setOfVariablesToReport):
    '''
    biomass (kg/m2)
    waterUseEfficiency biomass produced per volume of transpired water  kg/m3, typical 18
    maintenanceRate, maintenance (hour-1), typical 0.5 / (365*24)
    gamma, exponent relating runoff to maintenance for erosion
    alpha (half saturation constant relating grazing to vegetation density (kg/m2), typical 0.1
    dispersion (hour -1), dispersion rate, typical 0.2 / (365*24)
    sdOfNoise (kg/m2)
    timeStepDuration (hours)
    '''

    self.waterUseEfficiency=waterUseEfficiency
    self.biomass=biomass
    self.maintenanceRate=maintenanceRate
    self.runoff=runoff
    self.gamma=gamma
    self.alpha=alpha
    self.dispersion=dispersion
    self.sdOfNoise=sdOfNoise
    self.minimumAllowedBiomass=scalar(0.00000001)
    self.timeStepDuration=timeStepDuration
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport
    self.numberOfNeighbours=window4total(spatial(scalar(1)))

  def report(self,sample,timestep):
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                                'Xs': self.biomass,
                                'Xg': self.growth,
                                'Xm': self.maintenance,
                                'Xgr': self.grazing,
                                'Xdi': self.diffusion,
                                'Xno': self.noise,
                                'Xng': self.netGrowth
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = {
                                'Xs': self.biomass
                                }

    if timestep in self.timeStepsToReport:
      for variable in self.variablesToReport:
        report(self.variablesToReport[variable],generateNameST(variable, sample, timestep))


  def growthTerm(self,actualEvapotranspirationFlux):
    '''
    actualEvapotranspirationFlux, m/h, actual evapotranspiration
    uses waterUseEfficiency biomass produced per volume of transpired water  kg/m3, typical 18
    growth, kg m-2 h-1
    '''

    self.growth=self.waterUseEfficiency*actualEvapotranspirationFlux
    
  def maintenanceTerm(self):
    self.maintenance=scalar(0)-self.biomass*self.maintenanceRate

  def grazingTerm(self,grazingRate):
    '''
    grazingRate, kg m-2 h-1, typical 0.5 / (365*24)
    '''
    exponent=2.0
    self.grazing=scalar(0)-grazingRate*((self.biomass**exponent)/(self.biomass**exponent+self.alpha**exponent))

  def diffusionTerm(self):
    self.diffusion=self.dispersion*(window4total(spatial(scalar(self.biomass)))-self.numberOfNeighbours*self.biomass)

  def erosionTerm(self,runoff):
    '''
    runoff (overland flow per width) m3 h-1 m-1
    '''
    self.maintenanceErosion=runoff**self.gamma

  def noiseTerm(self):
    self.noise=normal(1)*self.sdOfNoise

  def netGrowthTerm(self,actualEvapotranspirationFlux,runoff,grazingRate):
    self.growthTerm(actualEvapotranspirationFlux)
    self.maintenanceTerm()
    self.grazingTerm(grazingRate)
    self.diffusionTerm()
    self.erosionTerm(runoff)
    self.noiseTerm()
    self.netGrowth=self.growth+self.maintenance+self.grazing+self.diffusion+self.maintenanceErosion+self.noise

  def update(self,actualEvapotranspirationFlux,runoff,grazingRate):
    self.netGrowthTerm(actualEvapotranspirationFlux,runoff,grazingRate)
    self.biomass=max(self.minimumAllowedBiomass,self.biomass+self.timeStepDuration*self.netGrowth)
    return self.biomass

