import math
from PCRaster import *
from PCRaster.Framework import *

# Notes:
# time step duration in h
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a waterslice over a whole cell)
# if unit cannot be derived in this way (e.g. flux/fluxAmount/store), unit is indicated
# inputs of function is PCRaster type, inside function Python types are used

# equations for the potential evapotranspiration Penman-Monteith
# based on the Thesis of van der Kwast, 2009

setclone('clone.map')

['airTemperature', 'albedo', 'clearSkyIncomingShortwaveRadiationFlatSurface', 'cloudFactor', 'elevationAboveSeaLevelOfMeteoStation', 'evapotranspirationOccurs', 'fWaterPotential', 'fractionReceivedFlatSurface', 'incomingShortwaveRadiation', 'incomingShortwaveRadiationFlatSurface', 'maxStomatalConduc', 'potentialEvapotranspirationAmount', 'potentialEvapotranspirationFlux', 'potentialEvapotranspirationFromCanopyAmount', 'potentialEvapotranspirationFromCanopyFlux', 'relativeHumidity', 'setOfVariablesToReport', 'timeStepDuration', 'timeStepsToReport', 'variablesToReport', 'vegHeight', 'windVelocityNotZero']

class EvapotranspirationPenman:
  def __init__(self, timeStepDuration,albedo,maxStomatalConductance, vegetationHeight, \
               timeStepsToReport,setOfVariablesToReport):

    # init only for suspend and resume in filter
    self.variablesToReport = {}
    self.airTemperature=scalar(0)
    self.clearSkyIncomingShortwaveRadiationFlatSurface=scalar(0)
    self.cloudFactor=scalar(0)
    self.elevationAboveSeaLevelOfMeteoStation=scalar(0)
    self.evapotranspirationOccurs=scalar(0)
    self.fWaterPotential=scalar(0)
    self.fractionReceivedFlatSurface=scalar(0)
    self.incomingShortwaveRadiation=scalar(0)
    self.incomingShortwaveRadiationFlatSurface=scalar(0)
    self.potentialEvapotranspirationAmount=scalar(0)
    self.potentialEvapotranspirationFlux=scalar(0)
    self.potentialEvapotranspirationFromCanopyAmount=scalar(0)
    self.potentialEvapotranspirationFromCanopyFlux=scalar(0)
    self.relativeHumidity=scalar(0)
    self.windVelocityNotZero=scalar(0)

    # real inits
    # maxStomatalConduc: maximum stomatal conductance (m s-1), appr. 0.04 m/s
    # vegetationHeight: height of vegetation, m
    self.timeStepDuration=scalar(timeStepDuration)
    self.albedo=scalar(albedo)
    self.maxStomatalConduc=maxStomatalConductance
    self.vegHeight=vegetationHeight

    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport

  def report(self,sample,timestep):
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                            'Ep': self.potentialEvapotranspirationFlux,
                            'Epc': self.potentialEvapotranspirationFromCanopyFlux,
                            'Ecl': self.cloudFactor,
                            'Ecs': self.clearSkyIncomingShortwaveRadiationFlatSurface,
                            'Eis': self.incomingShortwaveRadiationFlatSurface 
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = {}

    if timestep in self.timeStepsToReport:
      for variable in self.variablesToReport:
        report(self.variablesToReport[variable],generateNameST(variable, sample, timestep))


  def fluxToAmount(self,flux):
    fluxAmount=flux * self.timeStepDuration
    return fluxAmount

  def calculateClearSkySolarRadiation(self,fractionReceivedFlatSurface):
    # calculates incoming shortware radiation of a clear sky in W/m2, for a flat
    # surface (to make it comparable with actual incoming shortwave radiation measured at a meteo
    # station
    # solar constant (W/m2)taken from POTRAD manual and FAO Penman manual, Allen (no date), eg 28
    solarConstant=1368.0
    # equation 28, FAO Penman manual (Allen, no date), modified.
    extraTerrestrialRadiation=solarConstant*fractionReceivedFlatSurface
    # equation 37, FAO Penman manual (Allen, no date), modified.
    clearSkySolarRadiation=(0.75+0.00002*self.elevationAboveSeaLevelOfMeteoStation)*extraTerrestrialRadiation
    return clearSkySolarRadiation

  def potentialEvapotranspiration(self, \
                                  airTemperature, \
                                  relativeHumidity, \
                                  incomingShortwaveRadiation, \
                                  incomingShortwaveRadiationFlatSurface, \
                                  fractionReceivedFlatSurface, \
                                  windVelocity, \
                                  elevationAboveSeaLevelOfMeteoStation,
                                  fWaterPotential,
                                  evapotranspirationOccurs):
    # airTemperature (C), relativeHumidity (-), incomingShortRadiation (W m-2),
    # fractionReceivedFlatSurface (0-1),
    # incomingShortwaveRadiation, watt/m2 of soil surface (ie corrected for topography and solar angle, etc)
    # incomingShortwaveRadiationFlatSurface, watt/m2, of flat surface (ie raw data from meteo station)
    # fractionReceivedFlatSurface is the cross sectional area of a beam divided by its area of incidence on the soil surface
    # windVelocity, wind velocity, m/s
    # elevationAboveSeaLevelOfMeteoStation (m) is the elevation at the location where the
    # fWaterPotential (0-1) reduction factor for stomatal conductance
    self.airTemperature=scalar(airTemperature)
    self.relativeHumidity=scalar(relativeHumidity)
    self.incomingShortwaveRadiation=scalar(incomingShortwaveRadiation)
    self.incomingShortwaveRadiationFlatSurface=scalar(incomingShortwaveRadiationFlatSurface)
    self.fractionReceivedFlatSurface=scalar(fractionReceivedFlatSurface)
    self.windVelocityNotZero=max(0.000001,scalar(windVelocity))
    self.elevationAboveSeaLevelOfMeteoStation=scalar(elevationAboveSeaLevelOfMeteoStation)
    self.fWaterPotential=fWaterPotential
    self.evapotranspirationOccurs=evapotranspirationOccurs

    # constants
    ST=5.67*10.0**-8.0  # Stephan-Boltzman constant (W m-2 K-4)

    # eq 8.18 saturated vapour pressure (Pa), Allen et al 1998
    #sVapourPress=611.0*exp((17.27*self.airTemperature)/(237.3+self.airTemperature))
    sVapourPress=611.0*exp((17.27*self.airTemperature)/(237.3+self.airTemperature))

    # eq 8.17 actual vapour pressure (Pa)
    aVapourPress=self.relativeHumidity*sVapourPress

    self.clearSkyIncomingShortwaveRadiationFlatSurface=max(0.00000001,self.calculateClearSkySolarRadiation( \
                                                  self.fractionReceivedFlatSurface))

    # eq 8.19 cloud factor (-)
    self.cloudFactor=min(scalar(1.0), \
                     self.incomingShortwaveRadiationFlatSurface/self.clearSkyIncomingShortwaveRadiationFlatSurface)

    # eq 8.20 potential longwave radiation (W m-2), using Stefan-Boltzmann ST
    pLongRadiation=ST*(self.airTemperature+273.15)**4.0

    # eq 8.16 net long wave radiation (W m-2), Feddes et al 1983
    netLongRadiation=pLongRadiation*(0.56-0.008*sqrt(aVapourPress))*(0.1+0.9*self.cloudFactor)

    # eq 8.15 net radiation (W m-2)
    netRadiation=self.incomingShortwaveRadiation*(1.0-self.albedo)-netLongRadiation

    # aerodynamic resistance module
    # here with fixed input values for mesoWindHeight (m), meteoVegHeight (m), meteoWindHeight (m)

    # constants
    KAR=0.41          # von Karman constant

    # eq 8.22 displacement height (m), Rutter 1975
    displaHeight=0.75*self.vegHeight

    # eq 8.23 aerodynamic roughness (m), Rutter 1975
    aerodynamicRough=0.1*self.vegHeight

    # eq 8.24 wind extrapolator factor (-) for vegHeight>10
    mesoWindHeight=10.0 			# height of mesotropic wind (m)
    meteoVegHeight=0.1				# height of vegetation at meteo station (m)				
    meteoDisplaHeight=0.75*meteoVegHeight 	
    meteoAerodynamicRough=0.1*meteoVegHeight 	
    meteoWindHeight=3.0			# wind measurement height (m)
    windExtraFactor=ifthenelse(pcrlt(self.vegHeight,mesoWindHeight), 1.0,
         (ln((mesoWindHeight-meteoDisplaHeight)/meteoAerodynamicRough))/(ln((meteoWindHeight-meteoDisplaHeight)/meteoAerodynamicRough)))

    # eq 8.21 aerodynamic resistance (s m-1), using Von Karman constant KAR
    aerodynamicRes=ln(((mesoWindHeight-displaHeight)/aerodynamicRough)**2)/(KAR**2*self.windVelocityNotZero*windExtraFactor)

    # stomatal resistance module, simplified version (see Brolsma et al)

    stomatalConduc=max(self.maxStomatalConduc*self.fWaterPotential,0.00001)

    # stomatal resistance (s m-1)
    stomatalRes=1.0/stomatalConduc

    # other variables
    # here with fixed input values for Press (Pa)

    # contants
    airHeat=1013.0 		# air specific heat at constant pressure (J kg-1  K-1)
    latentHeat=2450000.0	# latent heat of water vaporization (J kg-1) from FAO56 Annex3

    # eq 8.29 vapour pressure deficit (Pa)
    defiVapourPress=sVapourPress-aVapourPress

    # eq 8.31 slope of saturation vapour pressure temperature relationship (Pa K-1) I didn't included 1/274.15 in Hans
    slopePressTemp=((4098.0*sVapourPress)/(237.3+self.airTemperature)**2.0)

    # eq 8.33 psychrometric constant (Pa K-1), FAO56 eq 8 (Brunt, 1956 see Annex3)
    Press=88400.0		# air pressure at z=1150 m (FAO56)
    psychrConstant=0.665*10.0**-3.0*Press

    # final equation for potential evapotranspiration

    # contants
    airDens=1.2047		# mean air density (kg m-3)

    # eq 8.14 potential evapotranspiration (mm s-1), Monteith, 1981 (or 1965?)
    pEvapt=(1.0/latentHeat)*(((netRadiation*slopePressTemp)+(airDens*airHeat*defiVapourPress/aerodynamicRes)) /
                  (slopePressTemp+psychrConstant*(1.0+stomatalRes/aerodynamicRes)))
    self.potentialEvapotranspirationFlux=pEvapt * (3600.0/1000.0)
    self.potentialEvapotranspirationFlux=ifthenelse( \
                                                   self.evapotranspirationOccurs, \
                                                   self.potentialEvapotranspirationFlux, \
                                                   scalar(0))
    # potential evapotranspiration amount (m over a timestep)
    self.potentialEvapotranspirationAmount=self.fluxToAmount(self.potentialEvapotranspirationFlux)

    # same with stomatal resistance zero (i.e., canopy evapotranspiration)
    stomatalRes=0.0
    pEvapt=(1.0/latentHeat)*(((netRadiation*slopePressTemp)+(airDens*airHeat*defiVapourPress/aerodynamicRes)) /
                  (slopePressTemp+psychrConstant*(1.0+stomatalRes/aerodynamicRes)))
    self.potentialEvapotranspirationFromCanopyFlux=pEvapt * (3600.0/1000.0)
    self.potentialEvapotranspirationFromCanopyFlux=ifthenelse( \
                                                   self.evapotranspirationOccurs, \
                                                   self.potentialEvapotranspirationFromCanopyFlux, \
                                                   scalar(0))
    # potential evapotranspiration amount (m over a timestep)
    self.potentialEvapotranspirationFromCanopyAmount=self.fluxToAmount(self.potentialEvapotranspirationFromCanopyFlux)

    return self.potentialEvapotranspirationFlux, self.potentialEvapotranspirationAmount,  \
           self.potentialEvapotranspirationFromCanopyFlux, self.potentialEvapotranspirationFromCanopyAmount

## test
#setclone('clone.map')
#timeStepDuration=scalar(1)
#albedo=scalar(0.3)
#maxStomatalConductance=0.0053
#vegetationHeight=8.0
#timeStepsToReportAll='test'
#setOfVariablesToReport='test'
#d_penmanPCRasterPython = EvapotranspirationPenman(
#                                         albedo, \
#                                         maxStomatalConductance, \
#                                         vegetationHeight, \
#                                         timeStepsToReportAll, \
#                                         setOfVariablesToReport)
## inputs for testing
#airTemperature=scalar(20)
#relativeHumidity=scalar(0.6)
#incomingShortwaveRadiationAtSurface=scalar(600.0)
#incomingShortwaveRadiationFlatSurface(600.0)
#fractionReceivedFlatSurface=1.0
#windVelocity=0.4
#elevationAboveSeaLevelOfMeteoStation=400.0
#fWaterPotential=0.2
## call function in class
#potentialEvapotranspirationFlux, potentialEvapotranspirationAmount, \
#potentialEvapotranspirationFromCanopyFlux, potentialEvapotranspirationFromCanopyAmount = \
#                            self.d_evapotranspirationPenman.potentialEvapotranspiration( \
#                            airTemperature, \
#                            relativeHumidity, \
#                            incomingShortwaveRadiationAtSurface, \
#                            incomingShortwaveRadiationFlatSurface, \
#                            fractionReceivedFlatSurface, \
#                            windVelocity,
#                            elevationAboveSeaLevelOfMeteoStation,
#                            fWaterPotential)
#report(potentialEvapotranspirationFlux,'testpetFlux.map')
#report(potentialEvapotranspirationAmount,'testpetAmount.map')
#report(potentialEvapotranspirationFromCanopyFlux,'testpetCanFlux.map')
#report(potentialEvapotranspirationFromCanopyAmount,'testpetCanAmount.map')
