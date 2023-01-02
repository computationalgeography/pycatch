from pcraster import *
import sys
import generalfunctions
from collections import deque
from pcraster.framework import *
import component

# notes
# time step duration in hours
# vertical fluxes, variable name 'flux'
#   water: m/hour (except where indicated)
#   vertical fluxes over a time step, variable name 'fluxAmount'
# amounts in storages, variable name 'store'
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

class BiomassModifiedMay(component.Component):
  def __init__(self, \
               biomass, \
               waterUseEfficiency, \
               maintenanceRate, \
               gamma, \
               alpha, \
               dispersion, \
               sdOfNoise, \
               LAIPerBiomass, \
               timeStepDuration,timeStepsToReport,setOfVariablesToReport):
    '''
    biomass (kg/m2)
        maximum values should be about
        1 kg/m2 (totaal, shrubland spain, prieto), maar lijken niet volgroeid
        2 kg/m2 (totaal, semi dry, south africa, abanda)
    waterUseEfficiency biomass produced per volume of transpired water  kg/m3, 
        18 (Siteur, MSc thesis)
        0.5 - 22 (Katerji) 
    maintenanceRate, maintenance (hour-1), typical 0.5 / (365*24)
        0.5 PER YEAR (Siteur, MSc thesis)
        verder niet makkelijk waarden te vinden
    gamma, exponent relating runoff to maintenance for erosion
        ?
    alpha (half saturation constant relating grazing to vegetation density (kg/m2), typical 0.4 or 0.5
    dispersion (hour -1), dispersion rate, typical 0.2 / (365*24)
    sdOfNoise (kg/m2 per hour??)
    LAIPerBiomass (/ kg/m2)
        LAI per biomass, e.g. use
        5.0 (LAI shrubs) / 2.0 (kg/m2, biomass of shrubs), is 2.5
    timeStepDuration (hours)
    '''

    self.waterUseEfficiency=waterUseEfficiency
    self.biomass=scalar(biomass)   # type definition seems to be required for reporting as maps
    self.maintenanceRate=maintenanceRate
    self.gamma=gamma
    self.alpha=alpha
    self.dispersion=dispersion
    self.sdOfNoise=sdOfNoise
    self.LAIPerBiomass=LAIPerBiomass
    self.LAI=scalar(0)
    self.growth=scalar(0)
    self.maintenance=scalar(0)
    self.grazing=scalar(0)
    self.diffusion=scalar(0)
    self.noise=scalar(0)
    self.netGrowth=scalar(0)
    #self.minimumAllowedBiomass=scalar(0.00000001)
    # let op higher minimum allowd biomass
    self.minimumAllowedBiomass=scalar(0.01)
    self.timeStepDuration=timeStepDuration
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport
    self.numberOfNeighbours=window4total(spatial(scalar(1)))

  def reportAsMaps(self, sample, timestep):
    self.output_mapping = {
                          'Xs': self.biomass,
                          'Xla': self.LAI,
                          'Xg': self.growth,
                          'Xm': self.maintenance,
                          'Xgr': self.grazing,
                          'Xdi': self.diffusion,
                          'Xno': self.noise,
                          'Xng': self.netGrowth
                           }
    self.variablesToReport = self.rasters_to_report(self.setOfVariablesToReport)
    self.reportMaps(sample, timestep)

  def growthTerm(self,actualEvapotranspirationFlux):
    '''
    actualEvapotranspirationFlux, m/h, actual evapotranspiration
    uses waterUseEfficiency biomass produced per volume of transpired water  kg/m3, typical 18
    growth, kg m-2 h-1
    '''

    self.growth=self.waterUseEfficiency*scalar(actualEvapotranspirationFlux)
    
  def maintenanceTerm(self):
    self.maintenance=scalar(0)-self.biomass*self.maintenanceRate

  def grazingTerm(self,grazingRate):
    '''
    grazingRate, kg m-2 h-1, typical 0.5 / (365*24)
    '''
    ## grazing met saturation
    #exponent=2.0 # gebruikte ik tot dec 2013, daarna dus exponent 1
    exponent=1.0
    self.grazing=scalar(0)-grazingRate*((self.biomass**exponent)/(self.biomass**exponent+self.alpha**exponent))

    #print 'other grazing lineair'
    #self.grazing=scalar(0)-grazingRate*self.biomass*2.0

    #print 'other grazing met hoek'
    #self.grazing=scalar(0) - max(0.0,0.0004*(self.biomass-grazingRate))
    #self.grazing=scalar(0) - max(0.0,0.00007*(self.biomass-grazingRate))

  def diffusionTerm(self):
    self.diffusion=self.dispersion*(window4total(spatial(scalar(self.biomass)))-self.numberOfNeighbours*self.biomass)

  def erosionTerm(self,runoff):
    '''
    runoff --  overland flow in m water depth per week 
               m water depth calculated as discharge in m3 divided by the cell area
    to get runoff m water depth per hour, runoff is divided by (7*24)
    so to get gamma using runoff in m water depth per hour, it has to be multiplied by (7*24)
    so for instance if gamma is 0.004, this would be 0.004*7*24=0.672
    if we would have 5 mm runoff in a week, this would be 0.005/(7*24)=2.9*10-5 m per hour
    this would give a maintenanceRateErosion of (with gamma is 0.004) 0.672*2.9*10-5=2*10-5
    note that this value is similar to the one Koen Siteur used for the whole maintenance rate, 5.7*10-5
    '''
    maintenanceRateErosion=self.gamma*runoff
    self.maintenanceErosion=scalar(0)-maintenanceRateErosion*self.biomass

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
    # for analysis only
    self.growthPart=self.growth+self.maintenance+self.diffusion+self.maintenanceErosion+self.noise

  def update(self,actualEvapotranspirationFlux,runoff,grazingRate):
    self.netGrowthTerm(actualEvapotranspirationFlux,runoff,grazingRate)
    self.biomass=max(self.minimumAllowedBiomass,self.biomass+self.timeStepDuration*self.netGrowth)
    self.LAI=self.LAIPerBiomass*self.biomass
    return self.biomass, self.LAI

  def setNewBiomass(self,biomass):
    self.biomass=biomass
    self.LAI=self.LAIPerBiomass*self.biomass


