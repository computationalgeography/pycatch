User manual
=============

Copyright of this software: Derek Karssenberg & Noemi Lana-Renault Monreal
Contact: d.karssenberg@uu.nl or noemi-solange.lana-renault@unirioja.es 


Installation
================

Create a conda environment for PCRaster (for instance by running pcraster_pycatch.yaml)

To run the model with 1 h timestep, run main.py
To run the model with 1 week timestep, run main_weekly.py (this model includes erosion, not documented yet)
To run the model with 10 second timestep, run main_seconds.py


Configuring the models
=============================

To remove all output, run clean.sh, this does NOT remove inputs, so there is no major risk, but check what is in clean.sh first

Settings are in configuration.py (1 h model) and configuration_weekly.py (1 week timestep model)


Instructions for hybrid modelling
====================================

Use the main_weekly.py model. Settings can be configured in configuration_weekly.py.

Realizations of the model are executed parallel (forked). Each realisation is written to a subdirectory. The number of
realizations can be defined in configuration_weekly.py. Each realisation of the model will be different as there is some tiny
random noise included which results in different paths of how the system evolves.

The model calculates grazing pressure inside the dynamic method. It currently gradually increases
grazing pressure then decreases it. Output is written to 1/grazing.npy, 2/grazing.npy,.. etc (all the same as there are no
differences between realizations.)

The model calculates and reports soil depth and biomass. These are written as txt files in the folders 1, 2, 3, etc again

The number of timesteps (time step duration is one week) is given in configuration_weekly.py

For displaying some map outputs, use setOfVariablesToReport = 'some' in configuration_weekly.py, it writes biomass and
grazing pressure for a defined interval. These can be displayed using the PCRaster Aguila visualisation tool.


Output file name conventions
----------------------------------

# first letter:

biomass or exchange vars starting with       X (Xst, biomass)
precipitation files starting with            P
interception files starting with             V (from vegetation)
surface store files starting with            S
infiltration files staring with              I
evapotranspiration files starting with       E
subsurface store files starting with         G (from groundwater)
runoff files starting with                   R
shading files starting with                  M (from shading due to Mountains )
budgets files starting with                  B
soil wash files starting with                W
regolith files starting with                 A (Ast, regolith depth (=soil depth))
bedrock weathering files starting with       C
base level files starting with               L
creep files starting with                    D

# second letter:

store (s, unit m)

actual flux in (i, unit m/h, for geomorphology m/year)
potential flux in (j, unit m/h, for geomorphology m/year)

actual flux out (o, unit m/h, for geomorphology m/year)
potential flux out (p, unit m/h, for geomorphology m/year)

actual flux in (in (positive) or out(negative)) (c, of change, unit m/h, for geomorphology m/year)
another flux (x, unit m/h, for geomorphology m/year)

lateral flux (q, cubic metre per hour)

two extra letters: other values, eg Ecl, cloud factor

see respective classes for details !!  E.g. actual flux out is not always
all out fluxes but sometimes only one of two.



For particle filtering (not yet tested)
-------------------------------------------

To create observed discharge, run createObservedDischarge.sh, it creates maps in the folder 'observations'. This is only required
for particle filtering

todo and changes are in changes.txt
at the bottom of changes is also the names of files used

file name conventions
~~~~~~~~~~~~~~~~~~~~~~

# reports as numpy arrays

Got.npy   self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour, from subsurface
Gxt.npy   self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour,from subsurface
Vot.npy   self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour, from canopy
Rq.npy    discharge m3/h
Rqs.npy   discharge m3/h, averaged over 2 hours

RPic.npy  self.maximumInterceptionCapacityPerLAI,
RPks.npy  self.ksat,
RPmm.npy  self.multiplierMaxStomatalConductance
RPrt.npy  self.regolithThicknessHomogeneous,
RPsc.npy  self.saturatedConductivityMetrePerDay,


Removed code
---------------

Early warning signals
~~~~~~~~~~~~~~~~~~~~~~

Code to calculate statistics in the dynamic section. Removed as it was hard to get GSTAT running.
Note also that it considerably slows down the model.

#    ############
#    # statistics
#    ############
#    # this needs GSTAT to run (rpy2) but cannot get it working
#    # SciKit Gstat (another tool may be useful as well) 
#    boundVector=(30.5,40.5)
#
#    # SOIL MOISTURE
#    self.d_subsurfaceWaterOneLayer.calculateSoilMoistureFraction()
#    variable=self.d_subsurfaceWaterOneLayer.soilMoistureFraction
#    variableSampled=ifthen(self.someLocs, variable)
#
#    self.historyOfSoilMoistureFraction=generalfunctions.keepHistoryOfMaps(self.historyOfSoilMoistureFraction, \
#                                                      variableSampled, \
#                                                      self.durationHistory)
#    stackOfMapsAsListVariable=list(self.historyOfSoilMoistureFraction)
#
#    if calculateStats:
#      if cfg.variances:
#        # temporal
#        dist,gamma=generalfunctions.experimentalVariogramValuesInTime(stackOfMapsAsListVariable,list(boundVector))
#        numpy.savetxt(generateNameST('sfT',self.currentSampleNumber(),self.currentTimeStep()),numpy.array(gamma))
#        # spatial
#        dist,gamma=generalfunctions.experimentalVariogramValues(stackOfMapsAsListVariable,boundVector,1,1,'test',2.0)
#        numpy.savetxt(generateNameST('sfS',self.currentSampleNumber(),self.currentTimeStep()),numpy.array(gamma))
#      # mean
#      #meanVariable=maptotal(variable)/self.numberOfCellsOnMap
#      meanVariable=areaaverage(variable,self.zoneMap)
#      generalfunctions.reportLocationsAsNumpyArray(self.aLocation,meanVariable,'sfA', \
#                       self.currentSampleNumber(),self.currentTimeStep())
#
#    # BIOMASS
#    variable=self.d_biomassModifiedMay.biomass
#    variableSampled=ifthen(self.someLocs, variable)
#
#    self.historyOfBiomass=generalfunctions.keepHistoryOfMaps(self.historyOfBiomass, \
#                                                      variableSampled, \
#                                                      self.durationHistory)
#    stackOfMapsAsListVariable=list(self.historyOfBiomass)
#
#    if calculateStats:
#      if cfg.variances:
#        # temporal
#        dist,gamma=generalfunctions.experimentalVariogramValuesInTime(stackOfMapsAsListVariable,list(boundVector))
#        numpy.savetxt(generateNameST('bioT',self.currentSampleNumber(),self.currentTimeStep()),numpy.array(gamma))
#        # spatial
#        dist,gamma=generalfunctions.experimentalVariogramValues(stackOfMapsAsListVariable,boundVector,1,1,'test',2.0)
#        numpy.savetxt(generateNameST('bioS',self.currentSampleNumber(),self.currentTimeStep()),numpy.array(gamma))
#      # mean
#      #meanVariable=maptotal(variable)/self.numberOfCellsOnMap
#      meanVariable=areaaverage(variable,self.zoneMap)
#      generalfunctions.reportLocationsAsNumpyArray(self.aLocation,meanVariable,'bioA', \
#                       self.currentSampleNumber(),self.currentTimeStep())
#
#    # REGOLITH THICKNESS
#    variable=self.d_regolithdemandbedrock.regolithThickness
#    variableSampled=ifthen(self.someLocs, variable)
#
#    self.historyOfRegolithThickness=generalfunctions.keepHistoryOfMaps(self.historyOfRegolithThickness, \
#                                                      variableSampled, \
#                                                      self.durationHistory)
#    stackOfMapsAsListVariable=list(self.historyOfRegolithThickness)
#
#    if calculateStats:
#      if cfg.variances:
#        # temporal
#        dist,gamma=generalfunctions.experimentalVariogramValuesInTime(stackOfMapsAsListVariable,list(boundVector))
#        numpy.savetxt(generateNameST('regT',self.currentSampleNumber(),self.currentTimeStep()),numpy.array(gamma))
#        # spatial
#        dist,gamma=generalfunctions.experimentalVariogramValues(stackOfMapsAsListVariable,boundVector,1,1,'test',2.0)
#        numpy.savetxt(generateNameST('regS',self.currentSampleNumber(),self.currentTimeStep()),numpy.array(gamma))
#      # mean
#      #meanVariable=maptotal(variable)/self.numberOfCellsOnMap
#      meanVariable=areaaverage(variable,self.zoneMap)
#      generalfunctions.reportLocationsAsNumpyArray(self.aLocation,meanVariable,'regA', \
#                       self.currentSampleNumber(),self.currentTimeStep())
#
#    # DEM 
#    variable=self.d_regolithdemandbedrock.dem
#    variableSampled=ifthen(self.someLocs, variable)
#
#    self.historyOfDem=generalfunctions.keepHistoryOfMaps(self.historyOfDem, \
#                                                      variableSampled, \
#                                                      self.durationHistory)
#    stackOfMapsAsListVariable=list(self.historyOfDem)
#
#    if calculateStats:
#      if cfg.variances:
#        # temporal
#        dist,gamma=generalfunctions.experimentalVariogramValuesInTime(stackOfMapsAsListVariable,list(boundVector))
#        numpy.savetxt(generateNameST('demT',self.currentSampleNumber(),self.currentTimeStep()),numpy.array(gamma))
#        # spatial
#        dist,gamma=generalfunctions.experimentalVariogramValues(stackOfMapsAsListVariable,boundVector,1,1,'test',2.0)
#        numpy.savetxt(generateNameST('demS',self.currentSampleNumber(),self.currentTimeStep()),numpy.array(gamma))
#      # mean
#      #meanVariable=maptotal(variable)/self.numberOfCellsOnMap
#      meanVariable=areaaverage(variable,self.zoneMap)
#      generalfunctions.reportLocationsAsNumpyArray( \
#                       self.aLocation,meanVariable,'demA',self.currentSampleNumber(),self.currentTimeStep())
#
#    # discharge 
#    downstreamEdge=generalfunctions.edge(self.clone,2,0)
#    pits=pcrne(pit(self.d_runoffAccuthreshold.ldd),0)
#    outflowPoints=pcrand(downstreamEdge,pits)
#    totQ=ifthen(self.clone,maptotal(ifthenelse(outflowPoints,self.d_runoffAccuthreshold.RunoffCubicMetrePerHour,scalar(0))))
#
#    variable=totQ
#    variableSampled=ifthen(self.someLocs, variable)
#
#    self.historyOfTotQ=generalfunctions.keepHistoryOfMaps(self.historyOfTotQ, \
#                                                      variableSampled, \
#                                                      self.durationHistory)
#    stackOfMapsAsListVariable=list(self.historyOfTotQ)
#
#    if calculateStats:
#      if cfg.variances:
#        # temporal
#        dist,gamma=generalfunctions.experimentalVariogramValuesInTime(stackOfMapsAsListVariable,list(boundVector))
#        numpy.savetxt(generateNameST('qT',self.currentSampleNumber(),self.currentTimeStep()),numpy.array(gamma))
#      # mean
#      generalfunctions.reportLocationsAsNumpyArray(self.aLocation,totQ,'qA',self.currentSampleNumber(),self.currentTimeStep())
#
#    # grazing rate
#    if calculateStats:
#      generalfunctions.reportLocationsAsNumpyArray( \
#                       self.aLocation,spatial(scalar(self.grazingRate)),'gA',self.currentSampleNumber(),self.currentTimeStep())

