import pcraster as pcr
import component

# Notes:
# time step duration in h
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a waterslice over a whole cell)
# if unit cannot be derived in this way (e.g. flux/fluxAmount/store), unit is indicated
# inputs of function is PCRaster type, inside function Python types are used

# equations for the potential evapotranspiration Penman-Monteith
# based on the Thesis of van der Kwast, 2009. Otherwise it is indicated.

['airTemperature', 'albedo', 'clearSkyIncomingShortwaveRadiationFlatSurface', 'cloudFactor', 'elevationAboveSeaLevelOfMeteoStation', 'evapotranspirationOccurs', 'fWaterPotential', 'fractionReceivedFlatSurface', 'incomingShortwaveRadiation', 'incomingShortwaveRadiationFlatSurface', 'maxStomatalConduc', 'potentialEvapotranspirationAmount', 'potentialEvapotranspirationFlux', 'potentialEvapotranspirationFromCanopyAmount', 'potentialEvapotranspirationFromCanopyFlux', 'relativeHumidity', 'setOfVariablesToReport', 'timeStepDuration', 'timeStepsToReport', 'variablesToReport', 'vegHeight', 'windVelocityNotZero']

class EvapotranspirationPenman(component.Component):
  def __init__(self, timeStepDuration, albedo, maxStomatalConductance, vegetationHeight, LAI,
               timeStepsToReport, setOfVariablesToReport):

    # init only for suspend and resume in filter
    self.variablesToReport = {}
    self.aboveVegHeight = pcr.scalar(0)
    self.aerodynamicRes = pcr.scalar(0)
    self.airTemperature = pcr.scalar(0)
    self.windVelocityAboveVegHeight = pcr.scalar(0)
    self.clearSkyIncomingShortwaveRadiationFlatSurface = pcr.scalar(0)
    self.cloudFactor = pcr.scalar(0)
    self.elevationAboveSeaLevelOfMeteoStation = pcr.scalar(0)
    self.evapotranspirationOccurs = pcr.scalar(0)
    self.fWaterPotential = pcr.scalar(0)
    self.fractionReceivedFlatSurface = pcr.scalar(0)
    self.incomingShortwaveRadiation = pcr.scalar(0)
    self.incomingShortwaveRadiationFlatSurface = pcr.scalar(0)
    self.potentialEvapotranspirationAmount = pcr.scalar(0)
    self.potentialEvapotranspirationFlux = pcr.scalar(0)
    self.potentialEvapotranspirationFromCanopyAmount = pcr.scalar(0)
    self.potentialEvapotranspirationFromCanopyFlux = pcr.scalar(0)
    self.relativeHumidity = pcr.scalar(0)
    self.windVelocityNotZero = pcr.scalar(0)

    # real inits
    # maxStomatalConduc: maximum stomatal conductance (m s-1), appr. 0.04 m/s
    # vegetationHeight: height of vegetation, m
    self.timeStepDuration = timeStepDuration
    self.albedo = pcr.scalar(albedo)
    self.maxStomatalConduc = maxStomatalConductance
    self.vegHeight = vegetationHeight
    self.LAI = LAI

    self.timeStepsToReport = timeStepsToReport
    self.setOfVariablesToReport = setOfVariablesToReport

  def reportAsMaps(self, sample, timestep):
    # reports
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                            'Ep': self.potentialEvapotranspirationFlux,
                            # 'gs': self.maxStomatalConduc,
                            # 'ra': self.aerodynamicRes,
                            # 'windz': self.windVelocityAboveVegHeight,
                            'Epc': self.potentialEvapotranspirationFromCanopyFlux,
                            # 'Ecl': self.cloudFactor,
                            # 'Ecs': self.clearSkyIncomingShortwaveRadiationFlatSurface,
                            # 'Eis': self.incomingShortwaveRadiationFlatSurface
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = { }
    self.reportMaps(sample, timestep)

  def updateVariablesAsNumpyToReport(self):
    self.variablesAsNumpyToReport = {
                            'gc': self.maxStomatalConduc,
                            # 'Ecl': self.cloudFactor,
                            # 'Ecs': self.clearSkyIncomingShortwaveRadiationFlatSurface,
                            # 'Eis': self.incomingShortwaveRadiationFlatSurface,
                                    }

  def reportAsNumpyOneFile(self, locations, sample, timestep, endTimeStep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpyOneFilePerRealization(locations, sample, timestep, endTimeStep)

  def reportAsNumpyMultipleFiles(self, locations, sample, timestep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpy(locations, sample, timestep)








  def fluxToAmount(self, flux):
    fluxAmount = flux * self.timeStepDuration
    return fluxAmount

  def calculateClearSkySolarRadiation(self, fractionReceivedFlatSurface):
    # calculates incoming shortware radiation of a clear sky in W/m2, for a flat
    # surface (to make it comparable with actual incoming shortwave radiation measured at a meteo
    # station
    # solar constant (W/m2)taken from POTRAD manual and FAO Penman manual, Allen (no date), eg 28
    solarConstant = 1368.0
    # equation 28, FAO Penman manual (Allen, no date), modified.
    extraTerrestrialRadiation = solarConstant * fractionReceivedFlatSurface
    # equation 37, FAO Penman manual (Allen, no date), modified.
    clearSkySolarRadiation = (0.75 + 0.00002 * self.elevationAboveSeaLevelOfMeteoStation) * extraTerrestrialRadiation
    return clearSkySolarRadiation

  def potentialEvapotranspiration(self,
                                  airTemperature,
                                  relativeHumidity,
                                  incomingShortwaveRadiation,
                                  incomingShortwaveRadiationFlatSurface,
                                  fractionReceivedFlatSurface,
                                  windVelocity,
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
    self.airTemperature = pcr.scalar(airTemperature)
    self.relativeHumidity = pcr.scalar(relativeHumidity)
    self.incomingShortwaveRadiation = pcr.scalar(incomingShortwaveRadiation)
    self.incomingShortwaveRadiationFlatSurface = pcr.scalar(incomingShortwaveRadiationFlatSurface)
    self.fractionReceivedFlatSurface = pcr.scalar(fractionReceivedFlatSurface)
    self.windVelocityNotZero = pcr.max(0.05, pcr.scalar(windVelocity)) #N FAO recommends u>0.5m/s
    self.elevationAboveSeaLevelOfMeteoStation = pcr.scalar(elevationAboveSeaLevelOfMeteoStation)
    self.fWaterPotential = fWaterPotential
    self.evapotranspirationOccurs = evapotranspirationOccurs


    # Radiation module

    ST = 5.67 * 10.0**-8.0  # Stephan-Boltzman constant (W m-2 K-4)

    # saturated vapour pressure (Pa), eq 8.18 (eq 11 in FAO56)
    sVapourPress = 611.0 * pcr.exp((17.27 * self.airTemperature) / (237.3 + self.airTemperature))

    # actual vapour pressure (Pa), eq 8.17
    aVapourPress = self.relativeHumidity * sVapourPress

    self.clearSkyIncomingShortwaveRadiationFlatSurface = pcr.max(0.00000001, self.calculateClearSkySolarRadiation(
                                                  self.fractionReceivedFlatSurface))

    # cloud factor (-), eq 8.19
    self.cloudFactor = pcr.min(pcr.scalar(1.0),
                     self.incomingShortwaveRadiationFlatSurface / self.clearSkyIncomingShortwaveRadiationFlatSurface)


    # potential longwave radiation (W m-2), eq 8.20
    pLongRadiation = ST * (self.airTemperature + 273.15)**4.0

    # net long wave radiation (W m-2) Feddes et al 1983, eq 8.16
    netLongRadiation = pLongRadiation * (0.56 - 0.008 * pcr.sqrt(aVapourPress)) * (0.1 + 0.9 * self.cloudFactor)

    # net radiation (W m-2), eq 8.15
    netRadiation = pcr.max(0.00001, self.incomingShortwaveRadiation * (1.0 - self.albedo) - netLongRadiation)


    # aerodynamic resistance module

    # Wind velocity at z (based on equation 5.33, Holton 2004)
    observedWindHeight = 3.0					# height of meteo station
    self.aboveVegHeight = self.vegHeight + 3.0				# we assume wind measurement is 3 m above vegetation
    self.windVelocityAboveVegHeight = self.windVelocityNotZero + (0.3 / 0.41) * pcr.ln(self.aboveVegHeight / observedWindHeight)


    #  aerodynamic resistance, Grace (1981)
    self.aerodynamicRes = pcr.ln((self.aboveVegHeight - 0.7 * self.vegHeight) / (0.1 * self.vegHeight))**2.0 / (0.41**2.0 * self.windVelocityAboveVegHeight)
    # self.aerodynamicRes=5.0
    # report(pcr.scalar(self.aerodynamicRes),'ra.map')


    # surface resistance module
    #  soil moisture reduction factor for stomatal conductance, simplified (see Brolsma et al 2010)

    stomatalConduc = pcr.max(self.maxStomatalConduc * self.fWaterPotential, 0.00001)

    # surface (or canopy) resistance (s m-1)
    stomatalRes = 1.0 / stomatalConduc
    LAIactive = 0.5 * self.LAI       # FAO Penman manual (Allen, 1998), Box 5

    surfaceRes = stomatalRes / LAIactive  # FAO Penman manual (Allen, 1998), equation 5

    # other variables
    # here with fixed input values for Press (Pa)

    # contants
    airHeat = 1013.0 		# air specific heat at constant pressure (J kg-1  K-1)
    latentHeat = 2450000.0	# latent heat of water vaporization (J kg-1) from FAO56 Annex3

    # eq 8.29 vapour pressure deficit (Pa)
    defiVapourPress = sVapourPress - aVapourPress

    # eq 8.31 slope of saturation vapour pressure temperature relationship (Pa K-1) I didn't included 1/274.15 in Hans
    slopePressTemp = ((4098.0 * sVapourPress) / (237.3 + self.airTemperature)**2.0)

    # eq 8.33 psychrometric constant (Pa K-1), FAO56 eq 8 (Brunt, 1956 see Annex3)
    Press = 88400.0		# air pressure at z=1150 m (FAO56)
    psychrConstant = 0.665 * 10.0**-3.0 * Press

    # final equation for potential evapotranspiration

    # contants
    airDens = 1.2047		# mean air density (kg m-3)

    # eq 8.14 potential evapotranspiration (mm s-1), Monteith, 1981 (or 1965?)
    pEvapt = (1.0 / latentHeat) * (((netRadiation * slopePressTemp) + (airDens * airHeat * defiVapourPress / self.aerodynamicRes)) /
                  (slopePressTemp + psychrConstant * (1.0 + surfaceRes / self.aerodynamicRes)))
    self.potentialEvapotranspirationFlux = pEvapt * (3600.0 / 1000.0)
    self.potentialEvapotranspirationFlux = pcr.ifthenelse(
                                                   self.evapotranspirationOccurs,
                                                   self.potentialEvapotranspirationFlux,
                                                   pcr.scalar(0))
    # potential evapotranspiration amount (m over a timestep)
    self.potentialEvapotranspirationAmount = self.fluxToAmount(self.potentialEvapotranspirationFlux)

    # same with surface resistance zero (i.e., canopy evapotranspiration)
    surfaceRes = 0.0
    pEvapt = (1.0 / latentHeat) * (((netRadiation * slopePressTemp) + (airDens * airHeat * defiVapourPress / self.aerodynamicRes)) /
                  (slopePressTemp + psychrConstant * (1.0 + surfaceRes / self.aerodynamicRes)))
    self.potentialEvapotranspirationFromCanopyFlux = pEvapt * (3600.0 / 1000.0)
    self.potentialEvapotranspirationFromCanopyFlux = pcr.ifthenelse(
                                                   self.evapotranspirationOccurs,
                                                   self.potentialEvapotranspirationFromCanopyFlux,
                                                   pcr.scalar(0))
    # potential evapotranspiration amount (m over a timestep)
    self.potentialEvapotranspirationFromCanopyAmount = self.fluxToAmount(self.potentialEvapotranspirationFromCanopyFlux)

    # report(self.potentialEvapotranspirationAmount,'jan')

    return self.potentialEvapotranspirationFlux, self.potentialEvapotranspirationAmount,  \
           self.potentialEvapotranspirationFromCanopyFlux, self.potentialEvapotranspirationFromCanopyAmount

# test
# setclone('mergeClone.map')
timeStepDuration = pcr.scalar(1)
albedo = pcr.scalar(0.15)
maxStomatalConductance = 0.0053
vegetationHeight = 5.0
timeStepsToReportAll = 'test'
setOfVariablesToReport = 'test'
# d_penmanPCRasterPython = EvapotranspirationPenman(
                                         # timeStepDuration, \
                                         # albedo, \
                                         # maxStomatalConductance, \
                                         # vegetationHeight, \
                                         # timeStepsToReportAll, \
                                         # setOfVariablesToReport)
# inputs for testing
airTemperature = pcr.scalar(23)
relativeHumidity = pcr.scalar(0.6)
incomingShortwaveRadiationAtSurface = pcr.scalar(300.0)
incomingShortwaveRadiationFlatSurface = pcr.scalar(300.0)
fractionReceivedFlatSurface = 1.0
windVelocity = 0.3
elevationAboveSeaLevelOfMeteoStation = 900.0
fWaterPotential = 0.2
evapotranspirationOccurs = 1


# call function in class
# potentialEvapotranspirationFlux, potentialEvapotranspirationAmount, \
# potentialEvapotranspirationFromCanopyFlux, potentialEvapotranspirationFromCanopyAmount = \
                            # d_penmanPCRasterPython.potentialEvapotranspiration( \
                            # airTemperature, \
                            # relativeHumidity, \
                            # incomingShortwaveRadiationAtSurface, \
                            # incomingShortwaveRadiationFlatSurface, \
                            # fractionReceivedFlatSurface, \
                            # windVelocity,
                            # elevationAboveSeaLevelOfMeteoStation,
                            # fWaterPotential, \
                            # evapotranspirationOccurs)

# report(potentialEvapotranspirationFlux,'testpetFlux.map')
# report(potentialEvapotranspirationAmount,'testpetAmount.map')
# report(potentialEvapotranspirationFromCanopyFlux,'testpetCanFlux.map')
# report(potentialEvapotranspirationFromCanopyAmount,'testpetCanAmount.map')



