from pcraster import *
import sys
from pcraster.framework import *
import component

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

class SoilWashMMF(component.Component):
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
    self.transportCapacityCParameter=((dFifty + 5.0)/0.32)**-0.6   # Hessel, Jetten, 2007 (from Govers)
    self.transportCapacityEtaParameter=((dFifty + 5.0)/300)**0.25
    self.criticalStreamPowerCmPerSec=0.4
    self.soilPorosity=soilPorosity
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport
    self.netDepositionKgPerCell=scalar(0)
    self.netDepositionMetre=scalar(0)
    self.transportCapacityVolumeFraction = scalar(0)
    self.streamPowerCmPerSec = scalar(0)
    self.accuDeposition=scalar(0)
    self.lateralFluxKg = scalar(0)
    self.totalDetachKgPerCell = scalar(0)
    self.transportCapacityKgPerCell = scalar(0)
    self.specificWeightRockKgPerCubicMetre=2650.0
    self.detachmentByRaindropsKgPerSquareMetre=scalar(0)
    self.flowVelocityMetrePerSecond=scalar(0)
    self.netDepositionMetreCum=scalar(0)

    # conversion kg rock to height soil
    volumeRockMaterialCubicMetresOfOneKiloRock=1.0/self.specificWeightRockKgPerCubicMetre
    volumeSoilOfOneKiloRock=volumeRockMaterialCubicMetresOfOneKiloRock/(1.0-self.soilPorosity)
    self.heightSoilOfOneKiloRock=volumeSoilOfOneKiloRock/cellarea()

  def reportAsMaps(self, sample, timestep):
    self.output_mapping = {
                           'Wde': self.netDepositionKgPerCell,
                           'Wdm': self.netDepositionMetre,
                           'Wfl': self.lateralFluxKg,
                           'Wdt': self.totalDetachKgPerCell,
                           'Wtc': self.transportCapacityKgPerCell,
                           'Wtv': self.transportCapacityVolumeFraction,
                           'Wsp': self.streamPowerCmPerSec,
                           'Wac': self.accuDeposition,
                           'Wdr': self.detachmentByRaindropsKgPerSquareMetre,
                           'Wdmc': self.netDepositionMetreCum
                          }
    self.variablesToReport = self.rasters_to_report(self.setOfVariablesToReport)
    self.reportMaps(sample, timestep)

  def updateStoneOrVegetationCover(self,vegetationCover):
    self.stoneOrVegetationCover=min(vegetationCover+self.stoneCoverFraction,1.0)

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
    #report(directThroughfallEnergy,'dtfe.map')
    return directThroughfallEnergy

  def kineticEnergyLeafDrainage(self):
    # LET OP: dit is totaal over een jaar, heb ik dus gedeeld door factor
    # deze component is behoorlijk klein tov kineticEnergy, wordt in VWK1 weggelaten
    leafDrainageEnergy=ifthenelse(self.plantHeightMetres < 0.15, scalar(0), (15.8 * sqrt(self.plantHeightMetres)) - 5.87)
    numberHoursRainPerYear=20.0*24.0
    numberHoursPerYear=365.0*24.0
    leafDrainageEnergyPerHour=leafDrainageEnergy*(numberHoursRainPerYear/numberHoursPerYear)
    leafDrainageEnergy=leafDrainageEnergyPerHour*self.rainstormDuration
    #report(leafDrainageEnergy,'lde.map')
    return leafDrainageEnergy

  def totalRainfallEnergy(self,erosiveRainfallIntensityFlux,directRainFlux):
    totalRainfallEnergy = self.kineticEnergyDirectRain(erosiveRainfallIntensityFlux,directRainFlux) + \
                          self.kineticEnergyLeafDrainage()
    return totalRainfallEnergy

  def detachmentSoilByRaindrops(self,erosiveRainfallIntensityFlux,directRainFlux):
    '''
    from Morgan & Duzant, 2008, Eqs 14-16
    assumes one grain size
    '''
    totalRainfallEnergy=self.totalRainfallEnergy(erosiveRainfallIntensityFlux,directRainFlux)
    self.detachmentByRaindropsKgPerSquareMetre=self.detachabilityOfSoilRaindrops * (1-self.stoneCoverFraction) * \
                                          totalRainfallEnergy * 0.001
    return self.detachmentByRaindropsKgPerSquareMetre


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
    '''
    This seems a combination of manning and A=alpha Q ^ Beta
    '''
    runoffCubicMetresPerSecond=runoffCubicMetresPerHour/(60.0*60.0)
    self.flowVelocityMetrePerSecond=( (runoffCubicMetresPerSecond**(2.0/5.0)) * (self.slope ** (3.0/10.0)) ) / \
                             ( (self.manningsN ** (3.0/5.0)) * (self.widthOfFlow ** (2.0 / 5.0)) )
    return self.flowVelocityMetrePerSecond
   
  def transportCapacity(self,runoffCubicMetresPerHour): 
    '''
    This is according to Hessel & Jetten 2007, the LISEM manual, and de Roo, 1996.
    They follow Govers, 1990.
    Conversion to cm is required, as the threshold in Govers 1990 is also in cm
    The unit of unit stream power is L/T, as it is calculated by multiplying
    flow velocity with gradient. So when it is calculated using velocity in m/s
    you get it in m/s and when it is calculated using velocity in cm/s, you get it in 
    cm/s.
    '''
    flowVelocityMetrePerSecond = self.flowVelocity(runoffCubicMetresPerHour)
    streamPowerMetrePerSecond = flowVelocityMetrePerSecond * self.slope
    self.streamPowerCmPerSec = (streamPowerMetrePerSecond * 100)
    self.transportCapacityVolumeFraction = self.transportCapacityCParameter * \
                                       max(0.0,(self.streamPowerCmPerSec - self.criticalStreamPowerCmPerSec)) ** \
                                       self.transportCapacityEtaParameter
    # original version up to feb 2014
    #transportCapacityKg=self.transportCapacityVolumeFraction*specificWeightKgPerCubicMetre*self.rainstormDuration
    # corrected version
    # transportCapacityVolumeFraction = S/(W+S) with S volume sand and W volume water, rewrite, giving S = FW/(1-F)
    self.transportCapacityCubicMetres=(self.transportCapacityVolumeFraction*(self.rainstormDuration*runoffCubicMetresPerHour))/ \
                                 (1.0-self.transportCapacityVolumeFraction)
    transportCapacityKg=self.transportCapacityCubicMetres*self.specificWeightRockKgPerCubicMetre
    return transportCapacityKg

  def kgToMetreHeight(self,kgSoilPerCell):
    heightSoil=self.heightSoilOfOneKiloRock*kgSoilPerCell
    return heightSoil 

  def calculateWash(self,runoffMetreWaterDepthPerHour,erosiveRainfallIntensityFlux,directRainFlux):
    '''
    directRainFlux           rain falling between gaps (i.e. through the gapfraction), represents the amount
                             i.e. it is a multiplier, it is the amount of water spreaded over the whole cell
    Note: in current main.py of pygeom this assumed to be net rain reaching soil surface, so water through gap fraction
    plus water that falls from leaves!!
    '''
    runoffDetachmentKgPerSquareMetre = self.detachmentSoilByRunoff(runoffMetreWaterDepthPerHour)
    raindropDetachmentKgPerSquareMetre = self.detachmentSoilByRaindrops(erosiveRainfallIntensityFlux,directRainFlux)

    self.totalDetachKgPerCell = cellarea()*(runoffDetachmentKgPerSquareMetre+raindropDetachmentKgPerSquareMetre)

    self.transportCapacityKgPerCell = self.transportCapacity(runoffMetreWaterDepthPerHour*cellarea())

    self.lateralFluxKg=accucapacityflux(self.ldd,self.totalDetachKgPerCell,self.transportCapacityKgPerCell)
    self.accuDeposition=accucapacitystate( \
                                self.ldd,self.totalDetachKgPerCell,self.transportCapacityKgPerCell)
    self.netDepositionKgPerCell=self.accuDeposition-self.totalDetachKgPerCell
    self.netDepositionMetre=self.kgToMetreHeight(self.netDepositionKgPerCell)
    self.netDepositionMetreCum=self.netDepositionMetre + self.netDepositionMetreCum
    return self.netDepositionKgPerCell, self.netDepositionMetre, self.lateralFluxKg, \
           self.totalDetachKgPerCell, self.transportCapacityKgPerCell

  def noWash(self):
    self.netDepositionKgPerCell=spatial(scalar(0))
    self.lateralFluxKg=spatial(scalar(0))
    self.totalDetachKgPerCell=spatial(scalar(0))
    self.transportCapacityKgPerCell=spatial(scalar(0))
    self.netDepositionMetre=spatial(scalar(0))
    return self.netDepositionKgPerCell, self.netDepositionMetre, self.lateralFluxKg, \
           self.totalDetachKgPerCell, self.transportCapacityKgPerCell

## test
#dem='../../switzerland/40m_small_area/mergeDem.map'
#ldd=lddcreate(dem,1e31,1e31,1e31,1e31)
#plantHeightMetres=5.0
#stoneCoverFraction=0.1
#vegetationCoverOfSoilFraction=0.1  
#manningsN=0.03 # 'original'
#detachabilityOfSoilRaindrops=1.6 # 'original'  (used for all scenarios)
#detachabilityOfSoilRunoff=6.4 #'original'
#durationOfRainstorm = 1.0
#soilPorosityFraction = 0.4
#   
#d_soilwashMMF=SoilWashMMF( \
#                                              ldd,
#                                              dem,
#                                              durationOfRainstorm,
#                                              plantHeightMetres,
#                                              detachabilityOfSoilRaindrops,
#                                              stoneCoverFraction, 
#                                              detachabilityOfSoilRunoff,
#                                              vegetationCoverOfSoilFraction,
#                                              manningsN,
#                                              soilPorosityFraction,
#                                              'tmp',
#                                              'tmp')
#
#rainfallFlux=0.005    # just the amount of rainfall given as a flux
#throughfallFlux=0.002
#runoffMetreWaterDepthPerHour=accuflux(ldd,rainfallFlux)
#
#netDeposition, netDepositionMetre, lateralFluxKg, totalDetachKgPerCell, transportCapacityKgPerCell= \
#                                            d_soilwashMMF.calculateWash( \
#                                            runoffMetreWaterDepthPerHour+0.000000001,rainfallFlux+0.000000001,throughfallFlux+0.000000001)
#
#report(netDeposition,'testnd')
#report(netDepositionMetre,'testndm')
#report(d_soilwashMMF.detachmentByRaindropsKgPerSquareMetre*cellarea(),'testdb') # raindrop detachment
#report(totalDetachKgPerCell,'testtdr')
#report(transportCapacityKgPerCell,'testtc')
#report(d_soilwashMMF.flowVelocityMetrePerSecond,'testfv')
