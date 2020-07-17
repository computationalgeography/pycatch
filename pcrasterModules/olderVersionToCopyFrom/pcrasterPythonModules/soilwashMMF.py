from PCRaster import *
import sys
from PCRaster.Framework import *

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

class SoilWashMMF:
  def __init__(self,ldd,dem,rainstormDuration,plantHeightMetres,detachabilityOfSoilRaindrops,stoneCoverFraction, \
               detachabilityOfSoilRunoff,vegetationCoverOfSoilFraction, manningsN,soilPorosity, \
               timeStepsToReport,setOfVariablesToReport):

    self.setSurfaceProperties(ldd,dem)

    self.rainstormDuration=rainstormDuration
    self.plantHeightMetres=plantHeightMetres
    self.detachabilityOfSoilRaindrops=detachabilityOfSoilRaindrops   # unit not clear J/m2 or g/J..
    self.stoneCoverFraction=stoneCoverFraction
    self.detachabilityOfSoilRunoff=detachabilityOfSoilRunoff
    self.vegetationCoverOfSoilFraction=vegetationCoverOfSoilFraction
    self.stoneOrVegetationCover=min(self.vegetationCoverOfSoilFraction+self.stoneCoverFraction,1.0)
    self.widthOfFlow=celllength()
    self.manningsN=manningsN
    dFifty=50.0  # mu metre median grain size
    self.transportCapacityCParameter=((dFifty + 5.0)/0.32)**-0.6
    self.transportCapacityEtaParameter=((dFifty + 5.0)/300)**0.25
    self.criticalStreamPowerCmPerSec=0.4
    self.soilPorosity=soilPorosity
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport
    self.netDepositionKgPerCell=scalar(0)
    self.netDepositionMetre=scalar(0)
    # conversion kg rock to height soil
    volumeRockMaterialCubicMetresOfOneKiloRock=1.0/2650.0
    volumeSoilOfOneKiloRock=volumeRockMaterialCubicMetresOfOneKiloRock/(1.0-self.soilPorosity)
    self.heightSoilOfOneKiloRock=volumeSoilOfOneKiloRock/cellarea()
  
  
  def report(self,sample,timestep):
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                                 'Wde': self.netDepositionKgPerCell,
                                 'Wdm': self.netDepositionMetre,
                                 'Wfl': self.lateralFluxKg,
                                 'Wdt': self.totalDetachKgPerCell,
                                 'Wtc': self.transportCapacityKgPerCell
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = {
                                 'Wde': self.netDepositionKgPerCell,
                                }

    if timestep in self.timeStepsToReport:
      for variable in self.variablesToReport:
        report(self.variablesToReport[variable],generateNameST(variable, sample, timestep))

  def setSurfaceProperties(self,ldd,dem):
    ''' needs to be updated when geom changes
    '''
    self.ldd=ldd
    self.dem=dem
    self.slope=slopeToDownstreamNeighbourNotFlat(self.dem,self.ldd,0.001)
    slopeAngle=atan(self.slope)
    self.sinOfSlope=sin(slopeAngle)


  def kineticEnergyDirectRain(self,erosiveRainfallIntensityFlux,directRainFlux):
    '''values of 8.95 and 8.44 are representative for UK
    erosiveRainfallIntensityFlux  preciptation flux, represents the intensity (power)
    directRainFlux                rain falling between gaps (i.e. through the gapfraction), represents the amount
                                  i.e. it is a multiplier, it is the amount of water spreaded over the whole cell
    '''
    
    erosiveRainfallIntensityMMPerHour=erosiveRainfallIntensityFlux*1000.0
    directRainAmountMM=directRainFlux*self.rainstormDuration*1000.0
    directThroughfallEnergyMayBecomeNegative=directRainAmountMM * \
                                             (8.95 + 8.44 * log10(erosiveRainfallIntensityMMPerHour))
    directThroughfallEnergy=max(0,directThroughfallEnergyMayBecomeNegative)
    report(directThroughfallEnergy,'dtfe.map')
    return directThroughfallEnergy

  def kineticEnergyLeafDrainage(self):
    # LET OP: dit is totaal over een jaar, heb ik dus gedeeld door factor
    # deze component is behoorlijk klein tov kineticEnergy, wordt in VWK1 weggelaten
    leafDrainageEnergy=ifthenelse(self.plantHeightMetres < 0.15, scalar(0), (15.8 * sqrt(self.plantHeightMetres)) - 5.87)
    numberHoursRainPerYear=20.0*24.0
    numberHoursPerYear=365.0*24.0
    leafDrainageEnergyPerHour=leafDrainageEnergy*(numberHoursRainPerYear/numberHoursPerYear)
    leafDrainageEnergy=leafDrainageEnergyPerHour*self.rainstormDuration
    report(leafDrainageEnergy,'lde.map')
    return leafDrainageEnergy

  def totalRainfallEnergy(self,erosiveRainfallIntensityFlux,directRainFlux):
    totalRainfallEnergy = self.kineticEnergyDirectRain(erosiveRainfallIntensityFlux,directRainFlux) + \
                          self.kineticEnergyLeafDrainage()
    return totalRainfallEnergy

  def detachmentSoilByRaindrops(self,erosiveRainfallIntensityFlux,directRainFlux):
    '''assume one grain size'''
    totalRainfallEnergy=self.totalRainfallEnergy(erosiveRainfallIntensityFlux,directRainFlux)
    detachmentByRaindropsKgPerSquareMetre=self.detachabilityOfSoilRaindrops * (1-self.stoneCoverFraction) * \
                                          totalRainfallEnergy * 0.001
    return detachmentByRaindropsKgPerSquareMetre


  def detachmentSoilByRunoff(self,runoffMetreWaterDepthPerHour):
    '''assume one grain size
    note that vegetation cover and stone cover of soil surface are summed 
    '''
    runoffMMWaterDepth=runoffMetreWaterDepthPerHour*1000*self.rainstormDuration
    detachmentByRunoffKgPerSquareMetre=self.detachabilityOfSoilRunoff * runoffMMWaterDepth**1.5 * \
                                       (1-self.stoneOrVegetationCover) * (self.sinOfSlope ** 0.3) \
                                       * 0.001
    return detachmentByRunoffKgPerSquareMetre

  def flowVelocity(self,runoffCubicMetresPerHour):
    runoffCubicMetresPerSecond=runoffCubicMetresPerHour/(60.0*60.0)
    flowVelocityMetrePerSecond=( (runoffCubicMetresPerSecond**(2.0/5.0)) * (self.slope ** (3.0/10.0)) ) / \
                             ( (self.manningsN ** (3.0/5.0)) * (self.widthOfFlow ** (2.0 / 5.0)) )
    return flowVelocityMetrePerSecond
   
  def transportCapacity(self,runoffCubicMetresPerHour): 
    flowVelocityMetrePerSecond = self.flowVelocity(runoffCubicMetresPerHour)
    streamPowerMetrePerSecond = flowVelocityMetrePerSecond * self.slope
    streamPowerCmPerSec = (streamPowerMetrePerSecond * 100)
    transportCapacityVolumeFraction = self.transportCapacityCParameter * \
                                       max(0.0,(streamPowerCmPerSec - self.criticalStreamPowerCmPerSec)) ** \
                                       self.transportCapacityEtaParameter
    specificWeightKgPerCubicMetre=2650
    transportCapacityKg=transportCapacityVolumeFraction*specificWeightKgPerCubicMetre*self.rainstormDuration
    return transportCapacityKg

  def kgToMetreHeight(self,kgSoilPerCell):
    heightSoil=self.heightSoilOfOneKiloRock*kgSoilPerCell
    return heightSoil 

  def calculateWash(self,runoffMetreWaterDepthPerHour,erosiveRainfallIntensityFlux,directRainFlux):
    runoffDetachmentKgPerSquareMetre = self.detachmentSoilByRunoff(runoffMetreWaterDepthPerHour)
    raindropDetachmentKgPerSquareMetre = self.detachmentSoilByRaindrops(erosiveRainfallIntensityFlux,directRainFlux)
    self.totalDetachKgPerCell = cellarea()*(runoffDetachmentKgPerSquareMetre+raindropDetachmentKgPerSquareMetre)
    self.transportCapacityKgPerCell = self.transportCapacity(runoffMetreWaterDepthPerHour*cellarea())
    self.lateralFluxKg=accucapacityflux(self.ldd,self.totalDetachKgPerCell,self.transportCapacityKgPerCell)
    self.netDepositionKgPerCell=accucapacitystate( \
                                self.ldd,self.totalDetachKgPerCell,self.transportCapacityKgPerCell)-self.totalDetachKgPerCell
    self.netDepositionMetre=self.kgToMetreHeight(self.netDepositionKgPerCell)
    return self.netDepositionKgPerCell, self.netDepositionMetre, self.lateralFluxKg, \
           self.totalDetachKgPerCell, self.transportCapacityKgPerCell

  def noWash(self):
    self.netDepositionKgPerCell=scalar(0)
    self.lateralFluxKg=scalar(0)
    self.totalDetachKgPerCell=scalar(0)
    self.transportCapacityKgPerCell=scalar(0)
    self.netDepositionMetre=scalar(0)
    return self.netDepositionKgPerCell, self.netDepositionMetre, self.lateralFluxKg, \
           self.totalDetachKgPerCell, self.transportCapacityKgPerCell

## test
#dem='mdtpaz4.map'
#ldd=lddcreate(dem,1e31,1e31,1e31,1e31)
#report(ldd,'ldd.map')
#rainstormDuration=1.0
#plantHeightMetres=5.0
#detachabilityOfSoilRaindrops=0.4
#detachabilityOfSoilRunoff=1.6
#stoneCoverFraction=0.1
#vegetationCoverOfSoilFraction=0.1
#erosiveRainfallIntensityFlux=0.001    # just the amount of rainfall given as a flux
#directRainFlux=0.0005          # 1-gapfraction value, should be smaller than erosiveRainfallIntensityFlux
#runoffMetreWaterDepthPerHour=accuflux(ldd,erosiveRainfallIntensityFlux/2.0)
#
#d_soilwashmmf = SoilWashMMF(ldd,dem,rainstormDuration,plantHeightMetres,detachabilityOfSoilRaindrops,stoneCoverFraction, \
#                            detachabilityOfSoilRunoff,vegetationCoverOfSoilFraction, 0.03)
#
##detachRunoffTonPerHectare=(detachRunoff*(100*100))/1000
##report(detachRunoff,'droton.map')
#
#netDepKg,latFluxKg,totDetachKg, transCapKg=d_soilwashmmf.calculateWash(runoffMetreWaterDepthPerHour,erosiveRainfallIntensityFlux,directRainFlux)
#
#report(netDepKg,'nd.map')
#report(latFluxKg,'lat.map')
#report(totDetachKg,'totdetach.map')
#report(transCapKg,'transcapkg.map')
#
#report(((netDepKg/2650)*(1.0/0.6))/cellarea(),'ndmetre.map')
