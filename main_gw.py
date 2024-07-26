# general
import datetime
import glob
import sys
import os

# add path to modules required
sys.path.append("./pcrasterModules/")

import configuration_gw as cfg

# PCRaster itself
import pcraster as pcr
import pcraster.framework as pcrfw


# from pcrasterModules
import datetimePCRasterPython
import interceptionuptomaxstore
import surfacestore
import infiltrationgreenandampt
import subsurfacewateronelayer_gw
import evapotranspirationpenman
import runoffaccuthreshold
import shading
import generalfunctions
import randomparameters

# from this folder
import exchangevariables


# only for advanced users
# uncomment following line and comment second line in case of particle filtering
# first time users should use class CatchmentModel(DynamicModel,MonteCarloModel):
# in case of particle filtering, CHANGE ALSO TIMESERIES FILE FOR SCENARIOS!!!!!!!!!
# class CatchmentModel(pcrfw.DynamicModel, pcrfw.MonteCarloModel, pcrfw.ParticleFilterModel):


class CatchmentModel(pcrfw.DynamicModel, pcrfw.MonteCarloModel):
  def __init__(self):
    pcrfw.DynamicModel.__init__(self)
    pcrfw.MonteCarloModel.__init__(self)
    pcr.setclone(cfg.cloneString)
    if cfg.filtering:
      pcrfw.ParticleFilterModel.__init__(self)

  def premcloop(self):
    self.clone = pcr.boolean(cfg.cloneString)
    self.dem = pcr.scalar(cfg.dem)
    self.createInstancesPremcloop()

    # required for reporting as numpy
    self.locations = pcr.cover(pcr.nominal(cfg.locations), 0)
    pcr.report(self.locations, 'locats')

    self.forestNoForest = pcr.boolean(cfg.forestNoForest)
    idMap = pcr.uniqueid(self.clone)
    oneLocationPerArea = pcr.areamaximum(idMap, pcr.spatial(self.forestNoForest)) == idMap
    self.locationsForParameters = pcr.cover(pcr.nominal(pcr.scalar(pcr.ifthen(oneLocationPerArea, self.forestNoForest)) + 1), 0)
    # end required for reporting as numpy

  def initial(self):
    self.timeStepDuration = cfg.timeStepDurationHoursFloatingPointValue
    self.initializeTime(cfg.startTimeYearValue, cfg.startTimeMonthValue, cfg.startTimeDayValue, self.timeStepDuration)
    self.createInstancesInitial()
    self.d_exchangevariables.upwardSeepageFlux = pcr.scalar(0)
    self.d_exchangevariables.evapFromSoilMultiplier = pcr.scalar(1)
    self.upwardSeepageFluxFromGroundWater = pcr.scalar(0)

    # budgets
    self.d_exchangevariables.cumulativePrecipitation = pcr.scalar(0)

  def dynamic(self):
    import generalfunctions

    # time
    self.d_dateTimePCRasterPython.update()
    timeDatetimeFormat = self.d_dateTimePCRasterPython.getTimeDatetimeFormat()

    # precipitation
    # for calibration
    rainfallFluxDeterm = pcr.timeinputscalar(cfg.rainfallFluxDetermTimeSeries, pcr.nominal(cfg.rainfallFluxDetermTimeSeriesAreas))
    # for the experiments
    rainfallFlux = rainfallFluxDeterm #generalfunctions.mapNormalRelativeError(rainfallFluxDeterm,0.25)
    self.d_exchangevariables.cumulativePrecipitation = \
            self.d_exchangevariables.cumulativePrecipitation + rainfallFlux * self.timeStepDuration

    # interception store
    actualAdditionFluxToInterceptionStore = self.d_interceptionuptomaxstore.addWater(rainfallFlux)
    throughfallFlux = rainfallFlux - actualAdditionFluxToInterceptionStore

    # surface store
    totalToSurfaceFlux = throughfallFlux + self.d_exchangevariables.upwardSeepageFlux + self.upwardSeepageFluxFromGroundWater
    potentialToSurfaceStoreFlux = self.d_surfaceStore.potentialToFlux()

    # potential infiltration
    potentialHortonianInfiltrationFlux = self.d_infiltrationgreenandampt.potentialInfiltrationFluxFunction()
    maximumSaturatedOverlandFlowInfiltrationFlux = self.d_subsurfaceWaterOneLayer.getMaximumAdditionFlux()
    potentialInfiltrationFlux = pcr.min(potentialHortonianInfiltrationFlux, maximumSaturatedOverlandFlowInfiltrationFlux)

    # abstraction from surface water
    potentialAbstractionFromSurfaceWaterFlux = potentialToSurfaceStoreFlux + potentialInfiltrationFlux
    actualAbstractionFromSurfaceWaterFlux, runoffCubicMetresPerHour = self.d_runoffAccuthreshold.update(
                                          totalToSurfaceFlux, potentialAbstractionFromSurfaceWaterFlux)
    potentialOutSurfaceStoreFlux = self.d_surfaceStore.potentialOutFlux()

    # infiltration
    availableForInfiltrationFlux = potentialOutSurfaceStoreFlux + actualAbstractionFromSurfaceWaterFlux
    availableForInfiltrationNotExceedingMaximumSaturatedOverlandFlowFlux = pcr.min(
                                          availableForInfiltrationFlux, maximumSaturatedOverlandFlowInfiltrationFlux)
    actualInfiltrationFlux = self.d_infiltrationgreenandampt.update(
                                          availableForInfiltrationNotExceedingMaximumSaturatedOverlandFlowFlux)

    # surface store
    surfaceStoreChange = actualAbstractionFromSurfaceWaterFlux - actualInfiltrationFlux
    self.d_surfaceStore.update(surfaceStoreChange)
    actualAdditionFlux = self.d_subsurfaceWaterOneLayer.addWater(actualInfiltrationFlux)

    if cfg.with_shading:
      # solar radiation (POTRAD, shading effect and inclination)
      fractionReceived, fractionReceivedFlatSurface, shaded = \
                                            self.d_shading.update(timeDatetimeFormat)

      # we assume all cells receive the same solar radiation as measured by the device
      # except for shading, if shading, there is nothing received
      fractionReceived = pcr.ifthenelse(shaded, pcr.scalar(0.0), pcr.scalar(1.0))
      #self.report(fractionReceived,'fr')
      #self.report(pcr.scalar(fractionReceivedFlatSurface),'frfs')
    else:
      fractionReceived = pcr.spatial(pcr.scalar(cfg.fractionReceivedValue))
      fractionReceivedFlatSurface = pcr.spatial(pcr.scalar(cfg.fractionReceivedFlatSurfaceValue))


    fWaterPotential = self.d_subsurfaceWaterOneLayer.getFWaterPotential()

    # potential evapotranspiration
    airTemperatureDeterm = pcr.timeinputscalar(cfg.airTemperatureDetermString, self.clone)
    airTemperature = airTemperatureDeterm #airTemperatureDeterm+mapnormal()

    relativeHumidityDeterm = pcr.timeinputscalar(cfg.relativeHumidityDetermString, self.clone)
    relativeHumidity = relativeHumidityDeterm #pcr.max(pcr.min(relativeHumidityDeterm+mapnormal()*0.1,pcr.scalar(1.0)),pcr.scalar(0))

    incomingShortwaveRadiationFlatSurface = pcr.timeinputscalar(cfg.incomingShortwaveRadiationFlatSurfaceString, self.clone)
    # incomingShortwaveRadiationFlatSurface = pcr.max(pcr.scalar(0),
    #                              generalfunctions.mapNormalRelativeError(incomingShortwaveRadiationFlatSurfaceDeterm,0.25))

    incomingShortwaveRadiationAtSurface = incomingShortwaveRadiationFlatSurface * fractionReceived

    windVelocityDeterm = pcr.timeinputscalar(cfg.windVelocityDetermString, self.clone)
    windVelocity = windVelocityDeterm #generalfunctions.mapNormalRelativeError(windVelocityDeterm,0.25)

    elevationAboveSeaLevelOfMeteoStation = cfg.elevationAboveSeaLevelOfMeteoStationValue

    potentialEvapotranspirationFlux, \
           potentialEvapotranspirationAmount, \
           potentialEvapotranspirationFromCanopyFlux, \
           potentialEvapotranspirationFromCanopyAmount = \
                            self.d_evapotranspirationPenman.potentialEvapotranspiration(
                            airTemperature,
                            relativeHumidity,
                            incomingShortwaveRadiationAtSurface,
                            incomingShortwaveRadiationFlatSurface,
                            fractionReceivedFlatSurface,
                            windVelocity,
                            elevationAboveSeaLevelOfMeteoStation,
                            fWaterPotential,
                            rainfallFlux < 0.000000000001)


    potentialEvapotranspirationFluxNoNegativeValues = pcr.max(0.0, potentialEvapotranspirationFlux)
    potentialEvapotranspirationFluxFromCanopyNoNegativeValues = pcr.max(0.0, potentialEvapotranspirationFromCanopyFlux)

    # evapotranspirate first from interception store
    actualAbstractionFluxFromInterceptionStore = self.d_interceptionuptomaxstore.abstractWater(
                            potentialEvapotranspirationFluxFromCanopyNoNegativeValues)

    # fraction of soil evapotranspiration depends on evapo from canopy
    evapFromSoilMultiplierMV = (potentialEvapotranspirationFluxFromCanopyNoNegativeValues -
                            actualAbstractionFluxFromInterceptionStore) / \
                            potentialEvapotranspirationFluxFromCanopyNoNegativeValues
    self.d_exchangevariables.evapFromSoilMultiplier = \
                           pcr.ifthenelse(potentialEvapotranspirationFluxNoNegativeValues < 0.0000000000001,
                           pcr.scalar(1), evapFromSoilMultiplierMV)

    # evapotranspirate from subsurface store
    # potentialEvapotranspirationFluxFromSubsurface= \
    #                       pcr.max(0.0,potentialEvapotranspirationFluxNoNegativeValues-actualAbstractionFluxFromInterceptionStore)
    potentialEvapotranspirationFluxFromSubsurface = self.d_exchangevariables.evapFromSoilMultiplier * \
                                                  potentialEvapotranspirationFluxNoNegativeValues
    actualAbstractionFluxFromSubsurface = self.d_subsurfaceWaterOneLayer.abstractWater(potentialEvapotranspirationFluxFromSubsurface)


    ## Ground water layer ##

    # potential percolation flux from soil layer, function of moisture content in soil layer only, it will not lead
    # to undersaturation of the soil layer, so potentialPercolation needs to check for this
    potentialPercolation = self.d_subsurfaceWaterOneLayer.potentialPercolation()

    # potential capillary rise from groundwater to soil layer, this is what the groundwater layer can give independent
    # of the soil layer status regarding storage status, again, it will not lead to undersaturation of the groundwater
    # layer, so potentialCapillaryRise needs to check for this
    unsaturatedConductivitySoilLayerMetrePerHour, saturationDegreeSoilLayer = \
                            self.d_subsurfaceWaterOneLayer.getUnsaturatedConductivity()
    potentialCapillaryRiseAmount = self.d_groundWaterLayer.getPotentialCapillaryRiseAmount( \
                             unsaturatedConductivitySoilLayerMetrePerHour, saturationDegreeSoilLayer)

    # percolation or capillary rise, depending on the soil moisture status of the two layers, water goes down
    # or up, not both
    saturationDegreeGroundWaterLayer = self.d_groundWaterLayer.updateDegreeOfSaturation()
    PercolationNoCapillaryRise = saturationDegreeSoilLayer > saturationDegreeGroundWaterLayer
    self.report(PercolationNoCapillaryRise, 'pncr')

    # take actual percolation from soil layer
    potentialPercolationIfPercolationOccurs = pcr.ifthenelse(PercolationNoCapillaryRise, potentialPercolation, pcr.scalar(0))
    actualPercolation = self.d_subsurfaceWaterOneLayer.abstractWater(potentialPercolationIfPercolationOccurs)

    # add actual percolation to the groundwater layer up to saturation, what is actually added is returned here 
    percolationFluxToGroundWater = self.d_groundWaterLayer.addWater(actualPercolation)

    # the groundwater may get saturated due to percolation and this water is considered upward seepage (routed down
    # over the land surface, note that this will be very small values
    upwardSeepageFluxFromPercolationToGroundwater = actualPercolation - percolationFluxToGroundWater

    # actual capillary rise
    potentialCapillaryRiseIfNoPercolationOccurs = pcr.ifthenelse( \
                              PercolationNoCapillaryRise, pcr.scalar(0), potentialCapillaryRiseAmount)
    actualCapillaryRise = self.d_subsurfaceWaterOneLayer.addWater(potentialCapillaryRiseIfNoPercolationOccurs)
    # input and output of function needs to be the same as groundwater conditions have already been checked for
    # above
    actualCapillaryRiseSecond = self.d_groundWaterLayer.abstractWater(actualCapillaryRise)
    self.report(actualPercolation-actualCapillaryRiseSecond,'vf')

    # upward seepage from groundwater
    # needs to be added to the rainfall/throughfall in the next timestep
    self.upwardSeepageFluxFromGroundWater = self.d_groundWaterLayer.lateralFlow() + \
                                       upwardSeepageFluxFromPercolationToGroundwater
    self.report(self.upwardSeepageFluxFromGroundWater,'usffg')

    ## End Ground water layer ##

    # upward seepage from subsurfacestore
    self.d_exchangevariables.upwardSeepageFlux = self.d_subsurfaceWaterOneLayer.lateralFlow()

    # reports
    self.reportComponentsDynamic()
    self.reportRandomParametersDynamic()
    self.printComponentsDynamic()
    if cfg.doReportComponentsDynamicAsNumpy:
      self.reportComponentsDynamicAsNumpy()

    #self.checkBudgets(self.currentSampleNumber(), self.currentTimeStep())

  def postmcloop(self):
    # required for reporting as numpy
    import generalfunctions
    self.timeStepDuration = cfg.timeStepDurationHoursFloatingPointValue # needed in case of forking, else the instances have been deleted
    self.initializeTime(cfg.startTimeYearValue, cfg.startTimeMonthValue, cfg.startTimeDayValue, self.timeStepDuration) # needed in case of forking, else the instances have been deleted
    self.createInstancesInitial() # needed in case of forking, else the instances have been deleted

    if cfg.doReportComponentsDynamicAsNumpy:
      self.reportAsNumpyComponentsPostmcloop()

  def createInstancesPremcloop(self):
    if cfg.with_shading:
      self.d_shading = shading.Shading(self.dem, cfg.latitudeOfCatchment, cfg.longitudeOfCatchment, cfg.timeZone, 1, cfg.timeStepsToReportRqs, cfg.shading_report_rasters)
      # print('no optimization of shading')


  def createInstancesInitial(self):
    import generalfunctions


    if cfg.readDistributionOfParametersFromDisk:
      path = '/home/derek/tmp/'
      maximumInterceptionCapacityPerLAI = pcr.scalar(path + pcrfw.generateNameS('RPic', self.currentSampleNumber()) + '.map')
      ksat = pcr.scalar(path + pcrfw.generateNameS('RPks', self.currentSampleNumber()) + '.map')
      regolithThicknessHomogeneous = pcr.scalar(path + pcrfw.generateNameS('RPrt', self.currentSampleNumber()) + '.map')
      saturatedConductivityMetrePerDay = pcr.scalar(path + pcrfw.generateNameS('RPsc', self.currentSampleNumber()) + '.map')
      multiplierMaxStomatalConductance = pcr.scalar(path + pcrfw.generateNameS('RPmm', self.currentSampleNumber()) + '.map')
    else:
      maximumInterceptionCapacityPerLAI = generalfunctions.areauniformBounds(
                                  0.0001, 0.0005, pcr.nominal(1), pcr.scalar(cfg.maximumInterceptionCapacityValue), cfg.createRealizations)
      ksat = generalfunctions.areauniformBounds(
                                  0.025, 0.05, pcr.nominal(1), pcr.scalar(cfg.ksatValue), cfg.createRealizations)
      regolithThicknessHomogeneous = generalfunctions.areauniformBounds(
                                  1.0, 3.5, cfg.areas, pcr.scalar(cfg.regolithThicknessHomogeneousValue), cfg.createRealizations)
      groundWaterLayerThicknessHomogeneous = generalfunctions.areauniformBounds(
                               1.0, 3.5, cfg.areas, pcr.scalar(cfg.groundWaterLayerThicknessHomogeneousValue), cfg.createRealizations)
      saturatedConductivityMetrePerDay = generalfunctions.mapuniformBounds(
                                  25.0, 40.0, pcr.scalar(cfg.saturatedConductivityMetrePerDayValue), cfg.createRealizations)
      multiplierMaxStomatalConductance = generalfunctions.mapuniformBounds(
                                  0.8, 1.1, pcr.scalar(cfg.multiplierMaxStomatalConductanceValue), cfg.createRealizations)

    if cfg.swapCatchments:
      regolithThicknessHomogeneous = generalfunctions.swapValuesOfTwoRegions(cfg.areas, regolithThicknessHomogeneous, True)

    self.d_randomparameters = randomparameters.RandomParameters(
                    cfg.timeStepsToReportRqs,
                    cfg.randomparameters_report_rasters,
                    maximumInterceptionCapacityPerLAI,
                    ksat,
                    regolithThicknessHomogeneous,
                    saturatedConductivityMetrePerDay,
                    multiplierMaxStomatalConductance)


    # class for exchange variables in initial and dynamic
    # introduced to make filtering possible
    self.d_exchangevariables = exchangevariables.ExchangeVariables(
                                    cfg.timeStepsToReportSome,
                                    cfg.exchange_report_rasters
                                    )

    ################
    # interception #
    ################

    self.ldd = cfg.lddMap

    initialInterceptionStore = pcr.scalar(0.000001)
    leafAreaIndex = pcr.scalar(cfg.leafAreaIndexValue)

    if cfg.swapCatchments:
      leafAreaIndex = generalfunctions.swapValuesOfTwoRegions(cfg.areas, leafAreaIndex, True)
    gapFraction = pcr.exp(-0.5 * leafAreaIndex)            # equation 40 in Brolsma et al 2010a
    maximumInterceptionStore = maximumInterceptionCapacityPerLAI * leafAreaIndex

    self.d_interceptionuptomaxstore = interceptionuptomaxstore.InterceptionUpToMaxStore(
                                    self.ldd,
                                    initialInterceptionStore,
                                    maximumInterceptionStore,
                                    gapFraction,
                                    cfg.calculateUpstreamTotals,
                                    self.timeStepDurationHours,
                                    cfg.timeStepsToReportSome,
                                    cfg.interception_report_rasters)

    #################
    # surface store #
    #################

    initialSurfaceStore = pcr.scalar(0.0)
    maxSurfaceStore = pcr.scalar(cfg.maxSurfaceStoreValue)
    self.d_surfaceStore = surfacestore.SurfaceStore(
                        initialSurfaceStore,
                        maxSurfaceStore,
                        self.timeStepDurationHours,
                        cfg.timeStepsToReportSome,
                        cfg.surfacestore_report_rasters)

    ################
    # infiltration #
    ################

    # N initialMoistureContentFraction taken from 1st July

    # DK
    # we do not use rts and Gs as input to calculate initial moisture fraction to avoid
    # problems when the initial regolith thickness is calibrated (it might be thinner than
    # initialMoistureThick -> problems!)
    # instead, we use initial moisture content fraction as input, read from disk, it is just calculated
    # by pcrcalc 'mergeInitialMoistureContentFraction=Gs000008.761/rts00008.761'
    # note that I also changed the name for the initial soil moisture as a fraction
    initialSoilMoistureFractionFromDisk = pcr.scalar(cfg.initialSoilMoistureFractionFromDiskValue)
    if cfg.swapCatchments:
      initialSoilMoistureFractionFromDisk = generalfunctions.swapValuesOfTwoRegions(cfg.areas, initialSoilMoistureFractionFromDisk, True)

    # initial soil moisture as a fraction should not be above soil porosity as a fraction, just a check
    soilPorosityFraction = pcr.scalar(cfg.soilPorosityFractionValue)
    if cfg.swapCatchments:
      soilPorosityFraction = generalfunctions.swapValuesOfTwoRegions(cfg.areas, soilPorosityFraction, True)
    initialSoilMoistureFraction = pcr.min(soilPorosityFraction, initialSoilMoistureFractionFromDisk)
    hf = pcr.scalar(-0.0000001)
    self.d_infiltrationgreenandampt = infiltrationgreenandampt.InfiltrationGreenAndAmpt(
                                    soilPorosityFraction,
                                    initialSoilMoistureFraction,
                                    ksat,
                                    hf,
                                    self.timeStepDurationHours,
                                    cfg.timeStepsToReportSome,
                                    cfg.infiltration_report_rasters)

    ####################
    # subsurface water #
    ####################

    demOfBedrockTopography = self.dem

    if cfg.mapsAsInput:
      stream = pcr.boolean(cfg.streamValue)
    else:
      upstreamArea = pcr.accuflux(cfg.lddMap,pcr.cellarea())
      stream = upstreamArea > 200000.0

    regolithThickness = pcr.ifthenelse(stream, 0.01, regolithThicknessHomogeneous)

    self.multiplierWiltingPoint = pcr.scalar(1.0)
    limitingPointFraction = pcr.scalar(cfg.limitingPointFractionValue)

    if cfg.swapCatchments:
      limitingPointFraction = generalfunctions.swapValuesOfTwoRegions(cfg.areas, limitingPointFraction, True)
    mergeWiltingPointFractionFS = pcr.scalar(cfg.mergeWiltingPointFractionFSValue)
    if cfg.swapCatchments:
      mergeWiltingPointFractionFS = generalfunctions.swapValuesOfTwoRegions(cfg.areas, mergeWiltingPointFractionFS, True)
    wiltingPointFractionNotChecked = mergeWiltingPointFractionFS * self.multiplierWiltingPoint
    wiltingPointFraction = pcr.min(wiltingPointFractionNotChecked, limitingPointFraction)

    fieldCapacityFraction = pcr.scalar(cfg.fieldCapacityFractionValue)
    if cfg.swapCatchments:
      fieldCapacityFraction = generalfunctions.swapValuesOfTwoRegions(cfg.areas, fieldCapacityFraction, True)

    self.d_subsurfaceWaterOneLayer = subsurfacewateronelayer_gw.SubsurfaceWaterOneLayer(
                                   self.ldd,
                                   demOfBedrockTopography,
                                   regolithThickness,
                                   initialSoilMoistureFraction,
                                   soilPorosityFraction,
                                   wiltingPointFraction,
                                   fieldCapacityFraction,
                                   limitingPointFraction,
                                   saturatedConductivityMetrePerDay,
                                   cfg.calculateUpstreamTotals,
                                   self.timeStepDurationHours,
                                   cfg.timeStepsToReportSome,
                                   cfg.subsurface_report_rasters,
                                   '')

    ####################
    # ground water #
    ####################

    # everything the same as subsurfacewateronelayer (soil) except:
    # thickness

    groundWaterLayerThickness = pcr.ifthenelse(stream, 0.2, groundWaterLayerThicknessHomogeneous)

    initialSoilMoistureFraction = fieldCapacityFraction


    self.d_groundWaterLayer = subsurfacewateronelayer_gw.SubsurfaceWaterOneLayer(
                                   self.ldd,
                                   demOfBedrockTopography,
                                   groundWaterLayerThickness,
                                   initialSoilMoistureFraction,
                                   soilPorosityFraction,
                                   wiltingPointFraction,
                                   fieldCapacityFraction,
                                   limitingPointFraction,
                                   saturatedConductivityMetrePerDay,
                                   cfg.calculateUpstreamTotals,
                                   self.timeStepDurationHours,
                                   cfg.timeStepsToReportSome,
                                   cfg.subsurface_report_rasters_gw,
                                   'G')


    ##########
    # runoff #
    ##########

    self.d_runoffAccuthreshold = runoffaccuthreshold.RunoffAccuthreshold(
                               self.ldd,
                               self.timeStepDurationHours,
                               cfg.timeStepsToReportRqs,
                               cfg.runoff_report_rasters)

    ######################
    # evapotranspiration #
    ######################

    albedo = pcr.scalar(cfg.albedoValue)
    if cfg.swapCatchments:
      albedo = generalfunctions.swapValuesOfTwoRegions(cfg.areas, albedo, True)

    maxStomatalConductance = pcr.scalar(cfg.maxStomatalConductanceValue) * multiplierMaxStomatalConductance
    if cfg.swapCatchments:
      maxStomatalConductance = generalfunctions.swapValuesOfTwoRegions(cfg.areas, maxStomatalConductance, True)

    vegetationHeight = pcr.scalar(cfg.vegetationHeightValue)
    if cfg.swapCatchments:
      vegetationHeight = generalfunctions.swapValuesOfTwoRegions(cfg.areas, vegetationHeight, True)
    self.d_evapotranspirationPenman = evapotranspirationpenman.EvapotranspirationPenman(
                                         self.timeStepDurationHours,
                                         albedo,
                                         maxStomatalConductance,
                                         vegetationHeight,
                                         leafAreaIndex,
                                         cfg.timeStepsToReportSome,
                                         cfg.evapotrans_report_rasters)


  def reportComponentsDynamic(self):
    """report dynamic components as PCRaster maps
    components, the modules that are reported
    see also reportAsNumpyComponentsPostmcloop
    """
    components = [ \
                  self.d_exchangevariables, \
                 self.d_randomparameters, \
                 self.d_interceptionuptomaxstore, \
                 # self.d_surfaceStore, \
                 self.d_infiltrationgreenandampt, \
                 self.d_evapotranspirationPenman, \
                 self.d_runoffAccuthreshold, \
                 self.d_subsurfaceWaterOneLayer, \
                 self.d_groundWaterLayer
                 ]

    if cfg.with_shading:
      components.append(self.d_shading)

    for component in components:
      component.reportAsMaps(self.currentSampleNumber(), self.currentTimeStep())


  def reportComponentsDynamicAsNumpy(self):
    """report dynamic components as PCRaster maps
    components, the modules that are reported
    see also reportAsNumpyComponentsPostmcloop
    """

    # report dynamic components as numpy, see also 'reportAsNupyComponentsPostmcloop'
    self.d_runoffAccuthreshold.reportAsNumpy(self.locations, self.currentSampleNumber(), self.currentTimeStep())
    self.d_subsurfaceWaterOneLayer.reportAsNumpy(self.locations, self.currentSampleNumber(), self.currentTimeStep())
    self.d_interceptionuptomaxstore.reportAsNumpy(self.locations, self.currentSampleNumber(), self.currentTimeStep())
    self.d_randomparameters.reportAsNumpy(self.locationsForParameters, self.currentSampleNumber(), self.currentTimeStep())


  def reportAsNumpyComponentsPostmcloop(self):
    """report dynamic components as PCRaster maps
    componentsToReportAsNumpy should correspond with the numpy one in reportComponentsDynamic
    """
    componentsToReportAsNumpy = [
                 self.d_runoffAccuthreshold,
                 self.d_subsurfaceWaterOneLayer,
                 self.d_interceptionuptomaxstore,
                 self.d_randomparameters
                                ]
    for component in componentsToReportAsNumpy:
      component.reportAsNumpyPostmcloop(range(1, cfg.nrOfSamples + 1), range(1, cfg.numberOfTimeSteps + 1))


  def reportRandomParametersDynamic(self):
    self.d_randomparameters.reportAtLastTimeStep(self.currentSampleNumber(), self.currentTimeStep(), self.nrTimeSteps())

  def printMemberVariables(self):
    import generalfunctions
    components = [
                 self.d_exchangevariables,
                 self.d_interceptionuptomaxstore,
                 self.d_surfaceStore,
                 self.d_infiltrationgreenandampt,
                 self.d_evapotranspirationPenman,
                 self.d_runoffAccuthreshold,
                 self.d_shading,
                 self.d_subsurfaceWaterOneLayer
                 ]

    for component in components:
      generalfunctions.printMemberVariables(component)

  def printComponentsDynamic(self):
    self.d_dateTimePCRasterPython.printit()

  def initializeTime(self, startTimeYear, startTimeMonth, startTimeDay, timeStepDurationHours):
    startTime = datetime.datetime(year=startTimeYear, month=startTimeMonth, day=startTimeDay)
    self.timeStepDurationHours = timeStepDurationHours
    self.timeStepDatetimeFormat = datetime.timedelta(hours=self.timeStepDurationHours)
    self.d_dateTimePCRasterPython = datetimePCRasterPython.DatetimePCRasterPython \
                                  (startTime, self.timeStepDatetimeFormat)

  def checkBudgets(self, currentSampleNumber, currentTimeStep):
    # DK not sure this is still correct

    increaseInPrecipitationStore = 0.0 - self.d_exchangevariables.cumulativePrecipitation
    pcr.report(increaseInPrecipitationStore, pcrfw.generateNameST('incP', currentSampleNumber, currentTimeStep))

    increaseInInterceptionStore = self.d_interceptionuptomaxstore.budgetCheck(currentSampleNumber, currentTimeStep)
    pcr.report(increaseInInterceptionStore, pcrfw.generateNameST('incI', currentSampleNumber, currentTimeStep))

    increaseInSurfaceStore = self.d_surfaceStore.budgetCheck(currentSampleNumber, currentTimeStep)
    pcr.report(increaseInSurfaceStore, pcrfw.generateNameST('incS', currentSampleNumber, currentTimeStep))
    increaseInSurfaceStoreQM = pcr.catchmenttotal(increaseInSurfaceStore, self.ldd) * pcr.cellarea()

    # let op: infiltration store is directly passed to subsurface store, thus is not a real store
    increaseInInfiltrationStore = self.d_infiltrationgreenandampt.budgetCheck(currentSampleNumber, currentTimeStep)

    increaseInSubSurfaceWaterStore, lateralFlowInSubsurfaceStore, abstractionFromSubSurfaceWaterStore = \
                             self.d_subsurfaceWaterOneLayer.budgetCheck(currentSampleNumber, currentTimeStep)
    increaseInSubSurfaceStoreQM = pcr.catchmenttotal(increaseInSubSurfaceWaterStore, self.ldd) * pcr.cellarea()

    increaseInRunoffStoreCubicMetresInUpstreamArea = self.d_runoffAccuthreshold.budgetCheck()

    totalIncreaseInStoresCubicMetresInUpstreamArea = 0.0
    stores = [increaseInPrecipitationStore, increaseInInterceptionStore, increaseInSurfaceStore, increaseInSubSurfaceWaterStore]
    for store in stores:
      increaseInStoreCubicMetresInUpstreamArea = pcr.catchmenttotal(store, self.ldd) * pcr.cellarea()
      totalIncreaseInStoresCubicMetresInUpstreamArea = totalIncreaseInStoresCubicMetresInUpstreamArea + \
                                                     increaseInStoreCubicMetresInUpstreamArea

    pcr.report(totalIncreaseInStoresCubicMetresInUpstreamArea, pcrfw.generateNameST('inSt', currentSampleNumber, currentTimeStep))
    pcr.report(increaseInRunoffStoreCubicMetresInUpstreamArea, pcrfw.generateNameST('inRu', currentSampleNumber, currentTimeStep))
    pcr.report(pcr.catchmenttotal(self.d_exchangevariables.upwardSeepageFlux, self.ldd) * pcr.cellarea(), pcrfw.generateNameST('inSe', currentSampleNumber, currentTimeStep))
    # total budget is total increase in stores plus the upward seepage flux for each ts that is passed to the next
    # timestep and thus not taken into account in the current timestep budgets
    budget = totalIncreaseInStoresCubicMetresInUpstreamArea + increaseInRunoffStoreCubicMetresInUpstreamArea + \
           lateralFlowInSubsurfaceStore * pcr.cellarea() + pcr.catchmenttotal(abstractionFromSubSurfaceWaterStore, self.ldd) * pcr.cellarea() + \
           pcr.catchmenttotal(self.d_exchangevariables.upwardSeepageFlux, self.ldd) * pcr.cellarea()
    pcr.report(budget, pcrfw.generateNameST('B-tot', currentSampleNumber, currentTimeStep))
    budgetRel = budget / increaseInRunoffStoreCubicMetresInUpstreamArea
    pcr.report(budgetRel, pcrfw.generateNameST('B-rel', currentSampleNumber, currentTimeStep))

  def suspend(self):
    import generalfunctions
    if self.currentTimeStep() != cfg.numberOfTimeSteps:
      self.timeStepForResume = self.currentTimeStep()

      components = [ self.d_exchangevariables,
                   self.d_randomparameters,
                   self.d_interceptionuptomaxstore,
                   self.d_surfaceStore,
                   self.d_infiltrationgreenandampt,
                   self.d_evapotranspirationPenman,
                   self.d_runoffAccuthreshold, \
                   # self.d_shading, \
                   self.d_subsurfaceWaterOneLayer]

      for component in components:
        generalfunctions.reportMemberVariablesOfAClassForSuspend(component, self.currentTimeStep(), self.currentSampleNumber())

  def updateWeight(self):
    print('#### UPDATEWEIGHTING')
    print('filter period', self.filterPeriod())
    print('filter timestep ', self._d_filterTimesteps[self.filterPeriod() - 1])
    print('lijst ', self._d_filterTimesteps)
    print('filter sample ', self.currentSampleNumber())
    modelledData = self.readmap('Rqs')
    observations = self.readDeterministic('observations/Rqs')
    #observations=pcr.ifthen(pit(self.ldd) != 0,syntheticData)
    measurementErrorSD = 3.0 * observations + 1.0
    sum = pcr.maptotal(((modelledData - observations)**2) / (2.0 * (measurementErrorSD**2)))
    weight = pcr.exp(-sum)
    weightFloatingPoint, valid = pcr.cellvalue(weight, 1)
    return weightFloatingPoint

  def resume(self):
    print('#### RESUMING')
    print('filter timesteps', self._d_filterTimesteps)
    print('filter period', self.filterPeriod())
    print('filter timestep', self._d_filterTimesteps[self.filterPeriod() - 2])

    import generalfunctions

    # rerun initial
    self.timeStepDuration = cfg.timeStepDurationHoursFloatingPointValue
    self.initializeTime(cfg.startTimeYearValue, cfg.startTimeMonthValue, cfg.startTimeDayValue, self.timeStepDuration)
    self.createInstancesInitial()
    self.d_exchangevariables.upwardSeepageFlux = pcr.scalar(0)
    self.d_exchangevariables.evapFromSoilMultiplier = pcr.scalar(1)
    self.d_exchangevariables.cumulativePrecipitation = pcr.scalar(0)

    # resume time information
    self.d_dateTimePCRasterPython.resume(self._d_filterTimesteps[self.filterPeriod() - 2])

    components = [ self.d_exchangevariables,
                  self.d_randomparameters,
                  self.d_interceptionuptomaxstore,
                  self.d_surfaceStore,
                  self.d_infiltrationgreenandampt,
                  self.d_evapotranspirationPenman,
                  self.d_runoffAccuthreshold, \
                  # self.d_shading, \
                  self.d_subsurfaceWaterOneLayer]

    for component in components:
      generalfunctions.readMemberVariablesOfAClassForResume(
                       component, self._d_filterTimesteps[self.filterPeriod() - 2], self.currentSampleNumber())

    # remove files used to resume
    for filename in glob.glob(str(self.currentSampleNumber()) + '/stateVar/*/*'):
      os.remove(filename)

if cfg.filtering:
  import generalfunctions
  myModel = CatchmentModel()
  dynamicModel = pcrfw.DynamicFramework(myModel, cfg.numberOfTimeSteps)
  mcModel = pcrfw.MonteCarloFramework(dynamicModel, cfg.nrOfSamples)
  mcModel.setForkSamples(True, 10)
  #pfModel = SequentialImportanceResamplingFramework(mcModel)
  pfModel = pcrfw.ResidualResamplingFramework(mcModel)
  filterTimestepsNoSelection = range(3750, cfg.numberOfTimeSteps + 1, 25)
  periodsToExclude = [
    [2617, 2976],
    [3649, 3689],
    [4173, 4416],
    [4046, 4366],
    [5281, 6075]
    ]
  filterTimesteps = generalfunctions.removePeriodsFromAListOfTimesteps(filterTimestepsNoSelection, periodsToExclude)
  pfModel.setFilterTimesteps(filterTimesteps)
  pfModel.run()

else:
  myModel = CatchmentModel()
  dynamicModel = pcrfw.DynamicFramework(myModel, cfg.numberOfTimeSteps)
  mcModel = pcrfw.MonteCarloFramework(dynamicModel, cfg.nrOfSamples)
  mcModel.setForkSamples(True, 10)
  mcModel.run()
