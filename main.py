# general
import math
import datetime
import random
import glob
import sys

# add path to modules required
sys.path.append("./pcrasterModules/")

# from pcrasterModules
import datetimePCRasterPython
import interceptionuptomaxstore
import surfacestore
import infiltrationgreenandampt
import subsurfacewateronelayer
import evapotranspirationpenman
import runoffaccuthreshold
import shading
import component
import generalfunctions

# from this folder
import exchangevariables
import randomparameters

# PCRaster itself
from pcraster import *
from pcraster.framework import *


# define number of hourly timesteps to run
numberOfTimeSteps = 10968
#numberOfTimeSteps = 1000

# Define the number of Monte Carlo samples or particles
# first time users will use 1 and results for that realization are written to
# the folder '1'
nrOfSamples = 1

# when classes of components are initialized, we pass a list with the time steps
# that are reported. These are defined here. In principle for each component
# a different set of time steps can be reported, by just passing another list
# but this version uses three different ones

# definition for components were all timesteps should be reported
timeStepsToReportAll = list(range(1,numberOfTimeSteps + 1,1))

# used for discharge only
timeStepsToReportRqs = list(range(20,numberOfTimeSteps + 1, 20))

# definition for components were a subset of timesteps should be reported
timeStepsToReportSome = list(range(100,numberOfTimeSteps + 1,100))

# switch to report for locations as small numpy files
# mainly used for particle filtering
doReportComponentsDynamicAsNumpy=False

# when True, a particle filtering run is done
# first time users should have this False
filtering=False

# selects whether a single, given, value is used for a number of parameters
# or whether a realization for that parameters is drawn
# first time users will use a single, fixed value for these parameters, so
# use False and search on createRealizations in the script to see which
# parameters are defined like this
createRealizations=False

# switch to swap parameter values between two catchments
# first time users will need to set this to False
swapCatchments=False

# when True, one can read a set of parameters for all Monte Carlo realizations
# from disk (e.g. representing probability distributions from a calibration)
# first time users should have a False here
readDistributionOfParametersFromDisk=False

# switch to define which set of variables are reported
# either 'full' or 'filtering'. These are passed to the class of a component
# where it the variables that are reported can be defined, i.e. either full or filtering
# early users will always use full
setOfVariablesToReport = 'full'
#setOfVariablesToReport = 'filtering'

# only for advanced users
# uncomment following line and comment second line in case of particle filtering
# first time users should use class CatchmentModel(DynamicModel,MonteCarloModel):
# in case of particle filtering, CHANGE ALSO TIMESERIES FILE FOR SCENARIOS!!!!!!!!!
#class CatchmentModel(DynamicModel,MonteCarloModel,ParticleFilterModel):
    

################
# model inputs #
################

# general ########

cloneString="mergeClone.map"
dem=scalar('mergeDem.map')
# report locations, i.e. outflow points, for instance, at the outlet
locations=nominal('mergeOutFlowsNominal.map')
# map with forest or no forest, only used when swapCatchments is True 
forestNoForest=boolean('mergeForestNoForest.map') 
areas='mergeForestNoForest.map'

# real time of first time step, duration of time step
# IMPORTANT NOTE: THIS IS NOW UTC TIME ALMOST CERTAINLY AT LEAST FOR SHADING
print("# IMPORTANT NOTE: THIS IS NOW UTC TIME ALMOST CERTAINLY AT LEAST FOR SHADING")
startTimeYearValue=2005
startTimeMonthValue=7
startTimeDayValue=1
timeStepDurationHoursFloatingPointValue=1.0  # only tested for one hour!!!!

# meteorology #######

# observed precipitation
rainfallFluxDetermTimeSeries="rainfallFluxTwoCatchsJulAugSep0506.tss"
# areas linked to rainfallFluxDetermTimeSeries
rainfallFluxDetermTimeSeriesAreas=nominal("mergeArnasSansaNominal.map")

airTemperatureDetermString="airTemperatureArnaJulAugSep0506.tss"
relativeHumidityDetermString="relativeHumidityArnasJulAugSep0506.tss"
incomingShortwaveRadiationFlatSurfaceString="incomingShortwaveRadiationArnasJulAugSep0506.tss"
windVelocityDetermString="windVelocityArnasJulAugSep0506.tss"
elevationAboveSeaLevelOfMeteoStationValue=900.0

# lat long for shading (solar radiation)
latitudeOfCatchment=52.12833333
longitudeOfCatchment=5.19861111
timeZone="Europe/Madrid"


# interception #######

maximumInterceptionCapacityValue=scalar(0.0002)
leafAreaIndexValue=scalar("mergeVegLAIFS.map")

# surface storage ######

maxSurfaceStoreValue=scalar(0.001)

# infiltration #######

# green and ampt
ksatValue=scalar(0.0163)
initialSoilMoistureFractionFromDiskValue=scalar("mergeFieldCapacityFractionFS.map")
soilPorosityFractionValue=scalar("mergePorosityFractionFS.map")


# regolith geometry ########

regolithThicknessHomogeneousValue=scalar(0.5)

# location of the stream, used to adjust regolith thickness there
streamValue=boolean('mergeStream.map')


# 'groundwater' (saturated flow) ##########

saturatedConductivityMetrePerDayValue=scalar(37.0)

limitingPointFractionValue=scalar("mergeLimitingPointFractionFS.map")
mergeWiltingPointFractionFSValue=scalar("mergeWiltingPointFractionFS.map")
fieldCapacityFractionValue=scalar("mergeFieldCapacityFractionFS.map")

# evapotranspiration ###########

# penman
multiplierMaxStomatalConductanceValue=scalar(1.0)
albedoValue=scalar("mergeVegAlbedoFS.map")
maxStomatalConductanceValue=scalar("mergeVegStomatalFS.map")
vegetationHeightValue=scalar("mergeVegHeightFS.map")



# dem geometry ###########

lddMap='mergeldd.map'





class CatchmentModel(DynamicModel,MonteCarloModel):
  def __init__(self):
    DynamicModel.__init__(self)
    MonteCarloModel.__init__(self)
    setclone(cloneString)
    if filtering:
      ParticleFilterModel.__init__(self)

  def premcloop(self):
    self.clone=boolean(cloneString)
    self.dem=dem
    self.createInstancesPremcloop()

    # required for reporting as numpy
    self.locations=cover(locations,0) 
    report(self.locations,'locats')
   
    self.forestNoForest=forestNoForest
    idMap=uniqueid(self.clone)
    oneLocationPerArea=pcreq(areamaximum(idMap,self.forestNoForest),idMap)
    self.locationsForParameters=cover(nominal(scalar(ifthen(oneLocationPerArea,self.forestNoForest))+1),0)
    # end required for reporting as numpy

  def initial(self):
    self.timeStepDuration = timeStepDurationHoursFloatingPointValue 
    self.initializeTime(startTimeYearValue,startTimeMonthValue,startTimeDayValue,self.timeStepDuration)
    self.createInstancesInitial()
    self.d_exchangevariables.upwardSeepageFlux=scalar(0)
    self.d_exchangevariables.evapFromSoilMultiplier=scalar(1)

    # budgets
    self.d_exchangevariables.cumulativePrecipitation=scalar(0)

  def dynamic(self):
    import generalfunctions

    # time
    self.d_dateTimePCRasterPython.update()
    timeDatetimeFormat=self.d_dateTimePCRasterPython.getTimeDatetimeFormat()
   
    # precipitation
    # for calibration
    rainfallFluxDeterm=timeinputscalar(rainfallFluxDetermTimeSeries, rainfallFluxDetermTimeSeriesAreas)
    ## for the experiments
    rainfallFlux = rainfallFluxDeterm #generalfunctions.mapNormalRelativeError(rainfallFluxDeterm,0.25)
    self.d_exchangevariables.cumulativePrecipitation= \
            self.d_exchangevariables.cumulativePrecipitation+rainfallFlux*self.timeStepDuration

    # interception store
    actualAdditionFluxToInterceptionStore=self.d_interceptionuptomaxstore.addWater(rainfallFlux)
    throughfallFlux=rainfallFlux-actualAdditionFluxToInterceptionStore

    # surface store
    totalToSurfaceFlux=throughfallFlux+self.d_exchangevariables.upwardSeepageFlux
    potentialToSurfaceStoreFlux=self.d_surfaceStore.potentialToFlux()

    # potential infiltration
    potentialHortonianInfiltrationFlux=self.d_infiltrationgreenandampt.potentialInfiltrationFluxFunction()
    maximumSaturatedOverlandFlowInfiltrationFlux=self.d_subsurfaceWaterOneLayer.getMaximumAdditionFlux()
    potentialInfiltrationFlux=min(potentialHortonianInfiltrationFlux,maximumSaturatedOverlandFlowInfiltrationFlux)

    # abstraction from surface water
    potentialAbstractionFromSurfaceWaterFlux=potentialToSurfaceStoreFlux + potentialInfiltrationFlux
    actualAbstractionFromSurfaceWaterFlux,runoffCubicMetresPerHour=self.d_runoffAccuthreshold.update( \
                                          totalToSurfaceFlux,potentialAbstractionFromSurfaceWaterFlux)
    potentialOutSurfaceStoreFlux=self.d_surfaceStore.potentialOutFlux()

    # infiltration
    availableForInfiltrationFlux=potentialOutSurfaceStoreFlux+actualAbstractionFromSurfaceWaterFlux
    availableForInfiltrationNotExceedingMaximumSaturatedOverlandFlowFlux=min( \
                                          availableForInfiltrationFlux,maximumSaturatedOverlandFlowInfiltrationFlux)
    actualInfiltrationFlux=self.d_infiltrationgreenandampt.update( \
                                          availableForInfiltrationNotExceedingMaximumSaturatedOverlandFlowFlux)

    # surface store
    surfaceStoreChange=actualAbstractionFromSurfaceWaterFlux-actualInfiltrationFlux
    self.d_surfaceStore.update(surfaceStoreChange)
    actualAdditionFlux=self.d_subsurfaceWaterOneLayer.addWater(actualInfiltrationFlux)

    # solar radiation (POTRAD, shading effect and inclination)
    fractionReceived, fractionReceivedFlatSurface,shaded= \
                                          self.d_shading.update(timeDatetimeFormat)

    # we assume all cells receive the same solar radiation as measured by the device
    # except for shading, if shading, there is nothing received
    fractionReceived=ifthenelse(shaded,scalar(0.0),scalar(1.0))
    
    fWaterPotential=self.d_subsurfaceWaterOneLayer.getFWaterPotential()

    # potential evapotranspiration
    airTemperatureDeterm=timeinputscalar(airTemperatureDetermString,self.clone)
    airTemperature = airTemperatureDeterm #airTemperatureDeterm+mapnormal() 

    relativeHumidityDeterm=timeinputscalar(relativeHumidityDetermString,self.clone)
    relativeHumidity = relativeHumidityDeterm #max(min(relativeHumidityDeterm+mapnormal()*0.1,scalar(1.0)),scalar(0))

    incomingShortwaveRadiationFlatSurface=timeinputscalar(incomingShortwaveRadiationFlatSurfaceString,self.clone)
    #incomingShortwaveRadiationFlatSurface = max(scalar(0),
    #                              generalfunctions.mapNormalRelativeError(incomingShortwaveRadiationFlatSurfaceDeterm,0.25))

    incomingShortwaveRadiationAtSurface=incomingShortwaveRadiationFlatSurface*fractionReceived

    windVelocityDeterm=timeinputscalar(windVelocityDetermString,self.clone)
    windVelocity = windVelocityDeterm #generalfunctions.mapNormalRelativeError(windVelocityDeterm,0.25)
    
    elevationAboveSeaLevelOfMeteoStation=elevationAboveSeaLevelOfMeteoStationValue
    
    potentialEvapotranspirationFlux, \
           potentialEvapotranspirationAmount, \
           potentialEvapotranspirationFromCanopyFlux, \
           potentialEvapotranspirationFromCanopyAmount = \
                            self.d_evapotranspirationPenman.potentialEvapotranspiration( \
                            airTemperature, \
                            relativeHumidity, \
                            incomingShortwaveRadiationAtSurface, \
                            incomingShortwaveRadiationFlatSurface, \
                            fractionReceivedFlatSurface, \
                            windVelocity, \
                            elevationAboveSeaLevelOfMeteoStation, \
                            fWaterPotential, \
                            pcrlt(rainfallFlux,0.000000000001))
 
    potentialEvapotranspirationFluxNoNegativeValues=max(0.0,potentialEvapotranspirationFlux)
    potentialEvapotranspirationFluxFromCanopyNoNegativeValues=max(0.0,potentialEvapotranspirationFromCanopyFlux)

    # evapotranspirate first from interception store
    actualAbstractionFluxFromInterceptionStore=self.d_interceptionuptomaxstore.abstractWater( \
                            potentialEvapotranspirationFluxFromCanopyNoNegativeValues)

    # fraction of soil evapotranspiration depends on evapo from canopy
    evapFromSoilMultiplierMV=(potentialEvapotranspirationFluxFromCanopyNoNegativeValues- \
                            actualAbstractionFluxFromInterceptionStore) / \
                            potentialEvapotranspirationFluxFromCanopyNoNegativeValues
    self.d_exchangevariables.evapFromSoilMultiplier= \
                           ifthenelse(potentialEvapotranspirationFluxNoNegativeValues < 0.0000000000001, \
                           scalar(1),evapFromSoilMultiplierMV)
                           
    # evapotranspirate from subsurface store
    #potentialEvapotranspirationFluxFromSubsurface= \
    #                       max(0.0,potentialEvapotranspirationFluxNoNegativeValues-actualAbstractionFluxFromInterceptionStore)
    potentialEvapotranspirationFluxFromSubsurface=self.d_exchangevariables.evapFromSoilMultiplier * \
                                                  potentialEvapotranspirationFluxNoNegativeValues
    actualAbstractionFluxFromSubsurface=self.d_subsurfaceWaterOneLayer.abstractWater(potentialEvapotranspirationFluxFromSubsurface)

    # upward seepage from subsurfacestore
    self.d_exchangevariables.upwardSeepageFlux=self.d_subsurfaceWaterOneLayer.lateralFlow()

    # reports
    self.reportComponentsDynamic()
    self.reportRandomParametersDynamic()
    self.printComponentsDynamic()
    if doReportComponentsDynamicAsNumpy:
      self.reportComponentsDynamicAsNumpy()

    #self.checkBudgets(self.currentSampleNumber(), self.currentTimeStep())

  def postmcloop(self):
    # required for reporting as numpy
    import generalfunctions
    self.timeStepDuration = timeStepDurationHoursFloatingPointValue # needed in case of forking, else the instances have been deleted
    self.initializeTime(startTimeYearValue,startTimeMonthValue,startTimeDayValue,self.timeStepDuration) # needed in case of forking, else the instances have been deleted
    self.createInstancesInitial() # needed in case of forking, else the instances have been deleted

    if doReportComponentsDynamicAsNumpy:
      self.reportAsNumpyComponentsPostmcloop()

  def createInstancesPremcloop(self):
    self.d_shading=shading.Shading(self.dem,latitudeOfCatchment,longitudeOfCatchment,timeZone,1,timeStepsToReportRqs,setOfVariablesToReport)
    #print 'no optimization of shading'
    

  def createInstancesInitial(self):
    import generalfunctions


    if readDistributionOfParametersFromDisk:
      path='/home/derek/tmp/'
      maximumInterceptionCapacityPerLAI=scalar(path + generateNameS('RPic', self.currentSampleNumber()) + '.map')
      ksat=scalar(path + generateNameS('RPks', self.currentSampleNumber()) + '.map')
      regolithThicknessHomogeneous=scalar(path + generateNameS('RPrt', self.currentSampleNumber()) + '.map')
      saturatedConductivityMetrePerDay=scalar(path + generateNameS('RPsc', self.currentSampleNumber()) + '.map')
      multiplierMaxStomatalConductance=scalar(path + generateNameS('RPmm', self.currentSampleNumber()) + '.map')
    else:
      maximumInterceptionCapacityPerLAI=generalfunctions.areauniformBounds( \
                                  0.0001,0.0005,nominal(1),maximumInterceptionCapacityValue,createRealizations)
      ksat=generalfunctions.areauniformBounds( \
                                  0.025,0.05,nominal(1),ksatValue,createRealizations)
      regolithThicknessHomogeneous=generalfunctions.areauniformBounds( \
                                  1.0,3.5,areas,regolithThicknessHomogeneousValue,createRealizations)
      saturatedConductivityMetrePerDay=generalfunctions.mapuniformBounds( \
                                  25.0,40.0,saturatedConductivityMetrePerDayValue,createRealizations)
      multiplierMaxStomatalConductance=generalfunctions.mapuniformBounds( \
                                  0.8,1.1,multiplierMaxStomatalConductanceValue,createRealizations)
       
    if swapCatchments:
      regolithThicknessHomogeneous=generalfunctions.swapValuesOfTwoRegions(areas,regolithThicknessHomogeneous,True)

    self.d_randomparameters=randomparameters.RandomParameters( \
                    timeStepsToReportRqs, \
                    setOfVariablesToReport, \
                    maximumInterceptionCapacityPerLAI, \
                    ksat,
                    regolithThicknessHomogeneous, \
                    saturatedConductivityMetrePerDay, \
                    multiplierMaxStomatalConductance)
  

    # class for exchange variables in initial and dynamic
    # introduced to make filtering possible
    self.d_exchangevariables=exchangevariables.ExchangeVariables( \
                                    timeStepsToReportSome, \
                                    setOfVariablesToReport, \
                                    )

    ################
    # interception #
    ################

    self.ldd=lddMap

    initialInterceptionStore=scalar(0.000001)
    leafAreaIndex=leafAreaIndexValue

    if swapCatchments:
      leafAreaIndex=generalfunctions.swapValuesOfTwoRegions(areas,leafAreaIndex,True)
    gapFraction=exp(-0.5*leafAreaIndex)            # equation 40 in Brolsma et al 2010a
    maximumInterceptionStore=maximumInterceptionCapacityPerLAI*leafAreaIndex

    self.d_interceptionuptomaxstore=interceptionuptomaxstore.InterceptionUpToMaxStore( \
                                    self.ldd, \
                                    initialInterceptionStore, \
                                    maximumInterceptionStore, \
                                    gapFraction, \
                                    self.timeStepDurationHours,
                                    timeStepsToReportSome,
                                    setOfVariablesToReport)

    #################
    # surface store #
    #################

    initialSurfaceStore=scalar(0.0)
    maxSurfaceStore=maxSurfaceStoreValue
    self.d_surfaceStore=surfacestore.SurfaceStore( \
                        initialSurfaceStore, \
                        maxSurfaceStore, \
                        self.timeStepDurationHours,
                        timeStepsToReportSome,
                        setOfVariablesToReport)

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
    initialSoilMoistureFractionFromDisk=initialSoilMoistureFractionFromDiskValue
    if swapCatchments:
      initialSoilMoistureFractionFromDisk=generalfunctions.swapValuesOfTwoRegions(areas,initialSoilMoistureFractionFromDisk,True)

    # initial soil moisture as a fraction should not be above soil porosity as a fraction, just a check
    soilPorosityFraction=soilPorosityFractionValue
    if swapCatchments:
      soilPorosityFraction=generalfunctions.swapValuesOfTwoRegions(areas,soilPorosityFraction,True)
    initialSoilMoistureFraction=min(soilPorosityFraction, initialSoilMoistureFractionFromDisk)
    hf=scalar(-0.0000001)
    self.d_infiltrationgreenandampt=infiltrationgreenandampt.InfiltrationGreenAndAmpt( \
                                    soilPorosityFraction, \
                                    initialSoilMoistureFraction, \
                                    ksat, \
                                    hf, \
                                    self.timeStepDurationHours,  \
                                    timeStepsToReportSome, \
                                    setOfVariablesToReport)

    ####################
    # subsurface water #
    ####################

    demOfBedrockTopography=self.dem

    stream=streamValue
    theSlope=slope(self.dem)
    regolithThickness=ifthenelse(stream,0.01,regolithThicknessHomogeneous)
  
    self.multiplierWiltingPoint=scalar(1.0)
    limitingPointFraction=limitingPointFractionValue

    if swapCatchments:
      limitingPointFraction=generalfunctions.swapValuesOfTwoRegions(areas,limitingPointFraction,True)
    mergeWiltingPointFractionFS=mergeWiltingPointFractionFSValue
    if swapCatchments:
      mergeWiltingPointFractionFS=generalfunctions.swapValuesOfTwoRegions(areas,mergeWiltingPointFractionFS,True)
    wiltingPointFractionNotChecked=mergeWiltingPointFractionFS*self.multiplierWiltingPoint
    wiltingPointFraction=min(wiltingPointFractionNotChecked,limitingPointFraction)

    fieldCapacityFraction=fieldCapacityFractionValue
    if swapCatchments:
      fieldCapacityFraction=generalfunctions.swapValuesOfTwoRegions(areas,fieldCapacityFraction,True)

    self.d_subsurfaceWaterOneLayer=subsurfacewateronelayer.SubsurfaceWaterOneLayer(
                                   self.ldd,
                                   demOfBedrockTopography,
                                   regolithThickness,
                                   initialSoilMoistureFraction,
                                   soilPorosityFraction,
                                   wiltingPointFraction,
                                   fieldCapacityFraction,
                                   limitingPointFraction,
                                   saturatedConductivityMetrePerDay,
                                   self.timeStepDurationHours,
                                   timeStepsToReportSome,
                                   setOfVariablesToReport)


    ##########
    # runoff #
    ##########

    self.d_runoffAccuthreshold=runoffaccuthreshold.RunoffAccuthreshold(
                               self.ldd,
                               self.timeStepDurationHours,
                               timeStepsToReportRqs,
                               setOfVariablesToReport)

    ######################
    # evapotranspiration #
    ######################

    albedo=albedoValue
    if swapCatchments:
      albedo=generalfunctions.swapValuesOfTwoRegions(areas,albedo,True)

    maxStomatalConductance=maxStomatalConductanceValue*multiplierMaxStomatalConductance 
    if swapCatchments:
      maxStomatalConductance=generalfunctions.swapValuesOfTwoRegions(areas,maxStomatalConductance,True)
    
    vegetationHeight=vegetationHeightValue
    if swapCatchments:
      vegetationHeight=generalfunctions.swapValuesOfTwoRegions(areas,vegetationHeight,True)
    self.d_evapotranspirationPenman=evapotranspirationpenman.EvapotranspirationPenman(
                                         self.timeStepDurationHours,\
                                         albedo, \
                                         maxStomatalConductance, \
                                         vegetationHeight, \
                                         leafAreaIndex, \
                                         timeStepsToReportSome, \
                                         setOfVariablesToReport)


  def reportComponentsDynamic(self):
    """report dynamic components as PCRaster maps
    components, the modules that are reported
    see also reportAsNumpyComponentsPostmcloop
    """
    components =[ \
                 #self.d_exchangevariables, \
                 self.d_randomparameters, \
                 self.d_interceptionuptomaxstore, \
                 #self.d_surfaceStore, \
                 self.d_infiltrationgreenandampt, \
                 self.d_evapotranspirationPenman, \
                 self.d_runoffAccuthreshold, \
                 self.d_shading, \
                 self.d_subsurfaceWaterOneLayer
                 ] 

    for component in components:
      component.reportAsMaps(self.currentSampleNumber(), self.currentTimeStep())


  def reportComponentsDynamicAsNumpy(self):
    """report dynamic components as PCRaster maps
    components, the modules that are reported
    see also reportAsNumpyComponentsPostmcloop
    """

    # report dynamic components as numpy, see also 'reportAsNupyComponentsPostmcloop'
    self.d_runoffAccuthreshold.reportAsNumpy(self.locations,self.currentSampleNumber(), self.currentTimeStep())
    self.d_subsurfaceWaterOneLayer.reportAsNumpy(self.locations,self.currentSampleNumber(), self.currentTimeStep())
    self.d_interceptionuptomaxstore.reportAsNumpy(self.locations,self.currentSampleNumber(), self.currentTimeStep())
    self.d_randomparameters.reportAsNumpy(self.locationsForParameters,self.currentSampleNumber(), self.currentTimeStep())


  def reportAsNumpyComponentsPostmcloop(self):
    """report dynamic components as PCRaster maps
    componentsToReportAsNumpy should correspond with the numpy one in reportComponentsDynamic
    """
    componentsToReportAsNumpy = [ \
                 self.d_runoffAccuthreshold, 
                 self.d_subsurfaceWaterOneLayer,
                 self.d_interceptionuptomaxstore,
                 self.d_randomparameters
                                ]
    for component in componentsToReportAsNumpy:
      component.reportAsNumpyPostmcloop(range(1,nrOfSamples+1),range(1,numberOfTimeSteps+1))
    

  def reportRandomParametersDynamic(self):
    self.d_randomparameters.reportAtLastTimeStep(self.currentSampleNumber(),self.currentTimeStep(),self.nrTimeSteps())

  def printMemberVariables(self):
    import generalfunctions
    components =[ \
                 self.d_exchangevariables, \
                 self.d_interceptionuptomaxstore, \
                 self.d_surfaceStore, \
                 self.d_infiltrationgreenandampt, \
                 self.d_evapotranspirationPenman, \
                 self.d_runoffAccuthreshold, \
                 self.d_shading, \
                 self.d_subsurfaceWaterOneLayer
                 ] 

    for component in components:
      generalfunctions.printMemberVariables(component)
    
  def printComponentsDynamic(self):
    self.d_dateTimePCRasterPython.printit()

  def initializeTime(self,startTimeYear, startTimeMonth, startTimeDay, timeStepDurationHours):
    startTime=datetime.datetime(year=startTimeYear, month=startTimeMonth, day=startTimeDay)
    self.timeStepDurationHours = timeStepDurationHours
    self.timeStepDatetimeFormat=datetime.timedelta(hours=self.timeStepDurationHours)
    self.d_dateTimePCRasterPython=datetimePCRasterPython.DatetimePCRasterPython \
                                  (startTime, self.timeStepDatetimeFormat)

  def checkBudgets(self,currentSampleNumber,currentTimeStep):

    increaseInPrecipitationStore=0.0-self.d_exchangevariables.cumulativePrecipitation
    report(increaseInPrecipitationStore,generateNameST('incP', currentSampleNumber, currentTimeStep))

    increaseInInterceptionStore=self.d_interceptionuptomaxstore.budgetCheck(currentSampleNumber, currentTimeStep)
    report(increaseInInterceptionStore,generateNameST('incI', currentSampleNumber, currentTimeStep))

    increaseInSurfaceStore=self.d_surfaceStore.budgetCheck(currentSampleNumber, currentTimeStep)
    report(increaseInSurfaceStore,generateNameST('incS', currentSampleNumber, currentTimeStep))
    increaseInSurfaceStoreQM=catchmenttotal(increaseInSurfaceStore,self.ldd)*cellarea()
    report(increaseInSurfaceStoreQM,generateNameST('testb', currentSampleNumber, currentTimeStep))

    # let op: infiltration store is directly passed to subsurface store, thus is not a real store
    increaseInInfiltrationStore=self.d_infiltrationgreenandampt.budgetCheck(currentSampleNumber, currentTimeStep)

    increaseInSubSurfaceWaterStore,lateralFlowInSubsurfaceStore,abstractionFromSubSurfaceWaterStore= \
                             self.d_subsurfaceWaterOneLayer.budgetCheck(currentSampleNumber , currentTimeStep)
    increaseInSubSurfaceStoreQM=catchmenttotal(increaseInSubSurfaceWaterStore,self.ldd)*cellarea()

    increaseInRunoffStoreCubicMetresInUpstreamArea=self.d_runoffAccuthreshold.budgetCheck()

    totalIncreaseInStoresCubicMetresInUpstreamArea=0.0
    stores=[increaseInPrecipitationStore, increaseInInterceptionStore, increaseInSurfaceStore, increaseInSubSurfaceWaterStore]
    for store in stores:
      increaseInStoreCubicMetresInUpstreamArea=catchmenttotal(store,self.ldd)*cellarea()
      totalIncreaseInStoresCubicMetresInUpstreamArea=totalIncreaseInStoresCubicMetresInUpstreamArea+ \
                                                     increaseInStoreCubicMetresInUpstreamArea
    
    report(totalIncreaseInStoresCubicMetresInUpstreamArea,generateNameST('inSt',currentSampleNumber,currentTimeStep))
    report(increaseInRunoffStoreCubicMetresInUpstreamArea ,generateNameST('inRu',currentSampleNumber,currentTimeStep))
    report(catchmenttotal(self.d_exchangevariables.upwardSeepageFlux,self.ldd)*cellarea(),generateNameST('inSe',currentSampleNumber,currentTimeStep))
    # total budget is total increase in stores plus the upward seepage flux for each ts that is passed to the next
    # timestep and thus not taken into account in the current timestep budgets
    budget=totalIncreaseInStoresCubicMetresInUpstreamArea+increaseInRunoffStoreCubicMetresInUpstreamArea + \
           lateralFlowInSubsurfaceStore*cellarea() + catchmenttotal(abstractionFromSubSurfaceWaterStore,self.ldd)*cellarea() + \
           catchmenttotal(self.d_exchangevariables.upwardSeepageFlux,self.ldd)*cellarea()
    report(budget,generateNameST('B-tot', currentSampleNumber, currentTimeStep))
    budgetRel=budget/increaseInRunoffStoreCubicMetresInUpstreamArea
    report(budgetRel,generateNameST('B-rel', currentSampleNumber, currentTimeStep))

  def suspend(self):
    import generalfunctions
    if self.currentTimeStep() != numberOfTimeSteps:
      self.timeStepForResume=self.currentTimeStep()

      components =[ self.d_exchangevariables, \
                   self.d_randomparameters, \
                   self.d_interceptionuptomaxstore, \
                   self.d_surfaceStore, \
                   self.d_infiltrationgreenandampt, \
                   self.d_evapotranspirationPenman, \
                   self.d_runoffAccuthreshold, \
                   #self.d_shading, \
                   self.d_subsurfaceWaterOneLayer] 

      for component in components:
        generalfunctions.reportMemberVariablesOfAClassForSuspend(component,self.currentTimeStep(),self.currentSampleNumber())

  def updateWeight(self):
    print('#### UPDATEWEIGHTING')
    print('filter period', self.filterPeriod())
    print('filter timestep ', self._d_filterTimesteps[self.filterPeriod()-1])
    print('lijst ', self._d_filterTimesteps)
    print('filter sample ', self.currentSampleNumber())
    modelledData=self.readmap('Rqs')
    observations=self.readDeterministic('observations/Rqs')
    #observations=ifthen(pit(self.ldd) != 0,syntheticData)
    measurementErrorSD=3.0*observations+1.0
    sum=maptotal(((modelledData-observations)**2)/(2.0*(measurementErrorSD**2)))
    weight=exp(-sum)
    weightFloatingPoint, valid=cellvalue(weight,1)
    return weightFloatingPoint

  def resume(self):
    print('#### RESUMING')
    print('filter timesteps', self._d_filterTimesteps)
    print('filter period', self.filterPeriod())
    print('filter timestep', self._d_filterTimesteps[self.filterPeriod()-2])

    import generalfunctions

    # rerun initial
    self.timeStepDuration = timeStepDurationHoursFloatingPointValue
    self.initializeTime(startTimeYearValue,startTimeMonthValue,startTimeDayValue,self.timeStepDuration)
    self.createInstancesInitial()
    self.d_exchangevariables.upwardSeepageFlux=scalar(0)
    self.d_exchangevariables.evapFromSoilMultiplier=scalar(1)
    self.d_exchangevariables.cumulativePrecipitation=scalar(0)

    # resume time information
    self.d_dateTimePCRasterPython.resume(self._d_filterTimesteps[self.filterPeriod()-2])

    components =[ self.d_exchangevariables, \
                  self.d_randomparameters, \
                  self.d_interceptionuptomaxstore, \
                  self.d_surfaceStore, \
                  self.d_infiltrationgreenandampt, \
                  self.d_evapotranspirationPenman, \
                  self.d_runoffAccuthreshold, \
                  #self.d_shading, \
                  self.d_subsurfaceWaterOneLayer] 

    for component in components:
      generalfunctions.readMemberVariablesOfAClassForResume( \
                       component,self._d_filterTimesteps[self.filterPeriod()-2],self.currentSampleNumber())

    # remove files used to resume
    for filename in glob.glob(str(self.currentSampleNumber()) + '/stateVar/*/*'):
      os.remove(filename)

if filtering:
  import generalfunctions
  myModel = CatchmentModel()
  dynamicModel = DynamicFramework(myModel,numberOfTimeSteps)
  mcModel = MonteCarloFramework(dynamicModel, nrOfSamples)
  mcModel.setForkSamples(True,10)
  #pfModel = SequentialImportanceResamplingFramework(mcModel)
  pfModel = ResidualResamplingFramework(mcModel)
  filterTimestepsNoSelection=range(3750,numberOfTimeSteps + 1,25)
  periodsToExclude=[ \
    [2617,2976],
    [3649,3689],
    [4173,4416],
    [4046,4366],
    [5281,6075]
    ]
  filterTimesteps=generalfunctions.removePeriodsFromAListOfTimesteps(filterTimestepsNoSelection,periodsToExclude)
  pfModel.setFilterTimesteps(filterTimesteps)
  pfModel.run()

else:
  myModel = CatchmentModel()
  dynamicModel = DynamicFramework(myModel, numberOfTimeSteps)
  mcModel = MonteCarloFramework(dynamicModel, nrOfSamples)
  mcModel.setForkSamples(True,10)
  mcModel.run()
