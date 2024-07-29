import pcraster as pcr
import pcraster.framework as pcrfw
import generalfunctions
import component

# Notes:
# time step duration in h
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# downward fluxes are positive
# amounts in storages in m (stores), variable name 'store'
# amounts as fraction (e.g. volumetric soil moisture content), variable name 'Fraction'
# amounts as waterlayer (i.e. volume fraction times thickness of layer, unit m), variable name 'Thick'
# (everything as a waterslice over a whole cell)
# if unit cannot be derived in this way (e.g. flux/fluxAmount/store), unit is indicated
# inputs of function is PCRaster type, inside function Python types are used

# the state of the block is given by the thicknesses, not the fractions, ie, thickness
# are always updated, fractions need to be retreived when needed

# geomorphology
# update regolith thickness and bedrock dem
# note that self.soilMoistureThick is the main state variable
# changing regolith thickness, update the following variables:

# conditions after changing regolith thickness:
# update self.soilPorosityThick after
# if self.soilMoistureThick < self.soilPorosityThick:
#   set self.soilMoistureThick to self.soilPorosityThick


# setclone('inputs/clone.map')

class SubsurfaceWaterOneLayer(component.Component):
  def __init__(self,
              ldd,
              demOfBedrockTopography,
              regolithThickness,
              initialSoilMoistureFraction,
              soilPorosityFraction,
              wiltingPointFraction,
              fieldCapacityFraction,
              limitingPointFraction,
              saturatedConductivityMetrePerDay,
              calculateUpstreamTotals,
              timeStepDuration,
              timeStepsToReport,
              setOfVariablesToReport,
              fileNamePrefix):

    ######
    # inits only for suspend and resume in filtering
    ######

    self.actualAbstractionFlux = pcr.scalar(0)
    self.actualAbstractionFluxAmount = pcr.scalar(0)
    self.actualAdditionFlux = pcr.scalar(0)
    self.actualAdditionFluxAmount = pcr.scalar(0)
    self.fWaterPotential = pcr.scalar(0)
    self.lateralFlowCubicMetrePerDay = pcr.scalar(0)
    self.lateralFlowCubicMetrePerTimeStep = pcr.scalar(0)
    self.lateralFlowDarcyFluxAmount = pcr.scalar(0)
    self.lateralFlowFluxAmount = pcr.scalar(0)
    self.lateralFlowFluxCubicMetresPerHour = pcr.scalar(0)
    self.maximumAbstractionThick = pcr.scalar(0)
    self.maximumAdditionThick = pcr.scalar(0)
    self.saturatedLayerThick = pcr.scalar(0)
    self.saturatedLayerThickness = pcr.scalar(0)
    self.soilMoistureMinusWiltingPointThick = pcr.scalar(0)
    self.upwardSeepageAmount = pcr.scalar(0)
    self.upwardSeepageFlux = pcr.scalar(0)
    self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour = pcr.scalar(0)
    self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour = pcr.scalar(0)
    self.totalSoilMoistureThickInUpstreamAreaCubicMetre = pcr.scalar(0)
    self.totalSoilMoistureFractionInUpstreamAreaTimesCellArea = pcr.scalar(0)
    self.totalSaturatedThickInUpstreamAreaCubicMetre = pcr.scalar(0)
    self.soilMoistureFraction = pcr.scalar(0)
    self.variablesToReport = {}
    self.variablesAsNumpyToReport = {}
    self.wiltingPointThick = pcr.scalar(0)
    self.potentialPercolationAmount = pcr.scalar(0)
    self.degreeOfSaturation = pcr.scalar(0)
    self.groundWaterDepthBelowSoil = pcr.scalar(0)
    self.potentialCapillaryRiseAmount = pcr.scalar(0)
    self.potentialCapillaryRiseNotCheckedFlux = pcr.scalar(0)
    self.unsaturatedConductivity = pcr.scalar(0)

    ####
    # real inits
    ####

    self.cellArea = pcr.cellarea()

    self.setBedrock(ldd, demOfBedrockTopography)

    # parameters depending on regolith thickness
    self.convertFractionsToThickness(regolithThickness, soilPorosityFraction, wiltingPointFraction,
                                  limitingPointFraction, fieldCapacityFraction)

    # main state variable, amount of water
    self.initialSoilMoistureFraction = initialSoilMoistureFraction
    self.initialSoilMoistureThick = self.initialSoilMoistureFraction * self.regolithThickness
    self.soilMoistureThick = self.initialSoilMoistureThick

    # values independent of regolith thickness
    self.saturatedConductivityMetrePerDay = saturatedConductivityMetrePerDay
    self.timeStepDuration = timeStepDuration
    self.timeStepsToReport = timeStepsToReport
    self.setOfVariablesToReport = setOfVariablesToReport
    self.calculateUpstreamTotals = calculateUpstreamTotals

    # DJ add check on fractions (e.g. wiltingpoint cannot be bigger than soil porosity)

    # for budgets
    self.actualAdditionCum = pcr.scalar(0)
    self.actualAbstractionCum = pcr.scalar(0)
    self.lateralFlowFluxAmountCum = pcr.scalar(0)
    self.upwardSeepageCum = pcr.scalar(0)

    # for filenames
    self.fileNamePrefix = fileNamePrefix

  def reportAsMaps(self, sample, timestep):
    self.output_mapping = {
                                self.fileNamePrefix + 'Gs': self.soilMoistureThick,
                                # 'Gad': self.maximumAdditionThick,
                                self.fileNamePrefix + 'Go': self.actualAbstractionFlux,
                                # 'Gi': self.actualAdditionFlux,
                                # 'Gq': self.lateralFlowFluxCubicMetresPerHour,
                                # 'Gwp': self.fWaterPotential,
                                # 'rts': self.regolithThickness,
                                # 'Gks': self.saturatedConductivityMetrePerDay
                                # 'Gxt': self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour,
                                # 'Got': self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour,
                                self.fileNamePrefix + 'Gppa': self.potentialPercolationAmount,
                                self.fileNamePrefix + 'Gdos': self.degreeOfSaturation,
                                self.fileNamePrefix + 'Ggwd': self.groundWaterDepthBelowSoil,
                                self.fileNamePrefix + 'Gpcr': self.potentialCapillaryRiseAmount,
                                self.fileNamePrefix + 'Gpcrf': self.potentialCapillaryRiseNotCheckedFlux,
                                self.fileNamePrefix + 'Gkun': self.unsaturatedConductivity
                           }
    # reports
    self.variablesToReport = self.rasters_to_report(self.setOfVariablesToReport)
    self.reportMaps(sample, timestep)

  def updateVariablesAsNumpyToReport(self):
    self.variablesAsNumpyToReport = {
                                'Gxt': self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour,
                                'Got': self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour,
                                'Gst': self.totalSoilMoistureThickInUpstreamAreaCubicMetre,
                                'Gtt': self.totalSoilMoistureFractionInUpstreamAreaTimesCellArea,
                                'GSt': self.totalSaturatedThickInUpstreamAreaCubicMetre
                                    }

  def reportAsNumpyOneFile(self, locations, sample, timestep, endTimeStep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpyOneFilePerRealization(locations, sample, timestep, endTimeStep)

  def reportAsNumpyMultipleFiles(self, locations, sample, timestep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpy(locations, sample, timestep)


  def convertFractionsToThickness(self, regolithThickness, soilPorosityFraction, wiltingPointFraction,
                                  limitingPointFraction, fieldCapacityFraction):
    # regolith thickness and parameters derived from regolith thickness,
    # need to be recalculated when regolith thickness changes
    self.regolithThickness = regolithThickness

    self.soilPorosityFraction = soilPorosityFraction
    self.soilPorosityThick = self.soilPorosityFraction * self.regolithThickness

    self.wiltingPointFraction = wiltingPointFraction
    self.wiltingPointThick = pcr.scalar(self.wiltingPointFraction * self.regolithThickness)

    self.limitingPointFraction = limitingPointFraction
    self.limitingPointThick = self.limitingPointFraction * self.regolithThickness

    self.fieldCapacityFraction = fieldCapacityFraction
    self.fieldCapacityThick = self.fieldCapacityFraction * self.regolithThickness

    self.maximumAvailableForVegetationUptakeThick = self.limitingPointThick - self.wiltingPointThick

  def setBedrock(self, ldd, demOfBedrockTopography):
    self.ldd = ldd
    self.demOfBedrockTopography = demOfBedrockTopography
    self.slopeToDownstreamNeighbourNotFlat = generalfunctions.slopeToDownstreamNeighbourNotFlat(
                                           self.demOfBedrockTopography, self.ldd, 0.001)


#  def report(self,sample,timestep):
#    self.variablesToReport = {}
#    if self.setOfVariablesToReport == 'full':
#      self.variablesToReport = {'Gx': self.upwardSeepageFlux,
#                                'Gs': self.soilMoistureThick,
#                                #'Gad': self.maximumAdditionThick,
#                                'Go': self.actualAbstractionFlux,
#                                #'Gi': self.actualAdditionFlux,
#                                #'Gq': self.lateralFlowFluxCubicMetresPerHour,
#                                #'Gwp': self.fWaterPotential,
#                                #'rts': self.regolithThickness,
#                                #'Gks': self.saturatedConductivityMetrePerDay
#                                'Gxt': self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour,
#                                'Got': self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour
#                                 }
#    if self.setOfVariablesToReport == 'filtering':
#      self.variablesToReport = {
#                                'Gsf': self.soilMoistureFraction
#                                }
#
#    if timestep in self.timeStepsToReport:
#      for variable in self.variablesToReport:
#        report(self.variablesToReport[variable],generateNameST(variable, sample, timestep))
#
#  def reportAsNumpy(self,locations,sample,timestep):
#    self.variablesAsNumpyToReport = {
#                                'Gxt': self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour,
#                                'Got': self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour,
#                                'Gst': self.totalSoilMoistureThickInUpstreamAreaCubicMetre,
#                                'Gtt': self.totalSoilMoistureFractionInUpstreamAreaTimesCellArea,
#                                'GSt': self.totalSaturatedThickInUpstreamAreaCubicMetre
#                                    }
#    import generalfunctions
#    for variable in self.variablesAsNumpyToReport:
#      generalfunctions.reportLocationsAsNumpyArray(locations,self.variablesAsNumpyToReport[variable], \
#                       variable,sample,timestep)
#
#  def reportAsNumpyPostmcloop(self,samples,timeSteps):
#    self.variablesAsNumpyToReport = {
#                                'Gxt': self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour,
#                                'Got': self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour,
#                                'Gst': self.totalSoilMoistureThickInUpstreamAreaCubicMetre,
#                                'Gtt': self.totalSoilMoistureFractionInUpstreamAreaTimesCellArea,
#                                'GSt': self.totalSaturatedThickInUpstreamAreaCubicMetre
#                                    }
#    import generalfunctions
#    for variable in self.variablesAsNumpyToReport:
#      aVariable = generalfunctions.openSamplesAndTimestepsAsNumpyArraysAsNumpyArray(variable,samples,timeSteps)
#      numpy.save(variable,aVariable)

  def fluxToAmount(self, flux):
    fluxAmount = flux * self.timeStepDuration
    return fluxAmount

  def amountToFlux(self, amount):
    flux = amount / self.timeStepDuration
    return flux

  def thicknessToFraction(self, thickness):
    fraction = thickness / self.regolithThickness
    return fraction

  def calculateMaximumAbstractionThick(self):
    # amount of water that can be taken out of the storage, up to wilting point
    self.maximumAbstractionThick = self.soilMoistureThick - self.wiltingPointThick

  def calculateMaximumAbstractionThickUpToFieldCapacity(self):
    # amount of water that can be taken out of the storage, up to wilting point
    self.maximumAbstractionThickUpToFieldCapacity = pcr.max(0.0, \
                                 self.soilMoistureThick - self.fieldCapacityThick)

  def getMaximumAbstractionFlux(self):
    self.calculateMaximumAbstractionThick()
    maximumAbstractionFlux = self.amountToFlux(self.maximumAbstractionThick)
    return maximumAbstractionFlux

  def calculateMaximumAdditionThick(self):
    # amount of water that can be added, up to saturation
    self.maximumAdditionThick = self.soilPorosityThick - self.soilMoistureThick

  def getMaximumAdditionFlux(self):
    self.calculateMaximumAdditionThick()
    maximumAdditionFlux = self.amountToFlux(self.maximumAdditionThick)
    return maximumAdditionFlux

  def abstractWater(self, potentialAbstractionFlux):
    # DJ add check on potentialAbstractionFlux >= 0
    self.calculateMaximumAbstractionThick()
    potentialAbstractionFluxAmount = self.fluxToAmount(potentialAbstractionFlux)
    self.actualAbstractionFluxAmount = pcr.min(self.maximumAbstractionThick, potentialAbstractionFluxAmount)
    self.soilMoistureThick = self.soilMoistureThick - self.actualAbstractionFluxAmount

    # conversions
    self.actualAbstractionFlux = self.amountToFlux(self.actualAbstractionFluxAmount)

    if self.calculateUpstreamTotals:
      self.totalActualAbstractionInUpstreamArea()

    return self.actualAbstractionFlux

  def addWater(self, potentialAdditionFlux):
    # DJ add check on potentialAdditionFlux >= 0
    self.calculateMaximumAdditionThick()
    potentialAdditionFluxAmount = self.fluxToAmount(potentialAdditionFlux)
    self.actualAdditionFluxAmount = pcr.min(self.maximumAdditionThick, potentialAdditionFluxAmount)
    self.soilMoistureThick = self.soilMoistureThick + self.actualAdditionFluxAmount

    # conversions
    self.actualAdditionFlux = self.amountToFlux(self.actualAdditionFluxAmount)

    return self.actualAdditionFlux

  def calculateThicknessOfSaturatedLayer(self):
    # thickness of saturated layer in m water depth (no pores)
    self.saturatedLayerThick = pcr.max(self.soilMoistureThick - self.fieldCapacityFraction * self.regolithThickness, pcr.scalar(0.0))
    # thickness of saturated layer in m (including pores)
    self.saturatedLayerThickness = self.saturatedLayerThick / self.soilPorosityFraction

  def calculateLateralFlowDarcy(self):
    self.calculateThicknessOfSaturatedLayer()
    self.lateralFlowCubicMetrePerDay = self.slopeToDownstreamNeighbourNotFlat * self.saturatedConductivityMetrePerDay \
                                * pcr.celllength() * self.saturatedLayerThickness;
    self.lateralFlowCubicMetrePerTimeStep = (self.lateralFlowCubicMetrePerDay / 24.0) * self.timeStepDuration
    self.lateralFlowDarcyFluxAmount = self.lateralFlowCubicMetrePerTimeStep / pcr.cellarea()

  def calculateSoilMoistureMinusWiltingPointThick(self):
    self.soilMoistureMinusWiltingPointThick = pcr.max(self.soilMoistureThick - self.wiltingPointThick, pcr.scalar(0.0))

  def calculateLateralFlow(self):
    self.calculateLateralFlowDarcy()
    # outgoing lateral flow may not be greater than soil moisture minus wilting point
    self.calculateSoilMoistureMinusWiltingPointThick()
    self.lateralFlowFluxAmount = pcr.min(self.lateralFlowDarcyFluxAmount, self.soilMoistureMinusWiltingPointThick)
    self.lateralFlowFluxCubicMetresPerHour = self.amountToFlux(self.lateralFlowFluxAmount) * pcr.cellarea()

  def lateralFlow(self):
    self.calculateLateralFlow()
    self.soilMoistureThick = self.soilMoistureThick - self.lateralFlowFluxAmount + pcr.upstream(self.ldd, self.lateralFlowFluxAmount)
    self.upwardSeepageAmount = pcr.max(pcr.scalar(0.0), self.soilMoistureThick - self.soilPorosityThick)
    self.soilMoistureThick = pcr.min(self.soilMoistureThick, self.soilPorosityThick)
    self.upwardSeepageFlux = self.amountToFlux(self.upwardSeepageAmount)

    if self.calculateUpstreamTotals:
      self.totalUpwardSeepageInUpstreamArea()
      self.totalSoilMoistureThickInUpstreamArea()
      self.totalSoilMoistureFractionInUpstreamArea()
      self.totalSaturatedThickInUpstreamArea()

    return self.upwardSeepageFlux

  def totalUpwardSeepageInUpstreamArea(self):
    self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour = pcr.accuflux(self.ldd, self.upwardSeepageFlux) * self.cellArea

  def totalActualAbstractionInUpstreamArea(self):
    self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour = pcr.accuflux(self.ldd, self.actualAbstractionFlux) * self.cellArea

  def totalSoilMoistureThickInUpstreamArea(self):
    self.totalSoilMoistureThickInUpstreamAreaCubicMetre = pcr.accuflux(self.ldd, self.soilMoistureThick) * self.cellArea

  def totalSoilMoistureFractionInUpstreamArea(self):
    self.calculateSoilMoistureFraction()
    self.totalSoilMoistureFractionInUpstreamAreaTimesCellArea = pcr.accuflux(self.ldd, self.soilMoistureFraction) * self.cellArea

  def totalSaturatedThickInUpstreamArea(self):
    self.totalSaturatedThickInUpstreamAreaCubicMetre = pcr.accuflux(self.ldd, self.saturatedLayerThick) * self.cellArea

  def getFWaterPotential(self):
    # calculates the reduction factor (0-1) for the stomatal conductance, i.e. stomatal
    # conductance in Penman will be the maximum stomatal conductance multiplied by
    # this reduction factor (see penman script)
    # the reduction factor is a linear function of soil moisture content, zero at
    # or below wilting point and one at or above field capacity, and a linear function inbetween
    fWaterPotentialTmp = (self.soilMoistureThick - self.wiltingPointThick) \
                           / self.maximumAvailableForVegetationUptakeThick
    fWaterPotentialTmpWilt = pcr.ifthenelse(self.soilMoistureThick < self.wiltingPointThick,
                           pcr.scalar(0), fWaterPotentialTmp)
    self.fWaterPotential = pcr.ifthenelse(self.soilMoistureThick > self.limitingPointThick,
                           pcr.scalar(1), fWaterPotentialTmpWilt)
    return self.fWaterPotential

  def updateDegreeOfSaturation(self):
    # degree of saturation calculated following van Beek 2009, used in
    # potentialPercolation method
    #degreeOfSaturationTmp = pcr.max(0.0,(self.soilMoistureThick-self.fieldCapacityThick)) / \
    #                        (self.soilPorosityThick-self.fieldCapacityThick)
    # idem, but using clapp 1979, which is somewhat different
    degreeOfSaturationTmp = self.soilMoistureThick/self.soilPorosityThick
    self.degreeOfSaturation = pcr.max(pcr.min(degreeOfSaturationTmp,pcr.scalar(1.0)), 0.0)
    return self.degreeOfSaturation

  def getUnsaturatedConductivity(self):
    # unsaturated conductivity of the layer, this is the same as the potential percolation
    # following Beek 2009, see also clapp 1979
    # we take b = 7, which is mean across textures, for sophistication use Clapp table 2
    # and adjust according to texture
    tmp  = self.updateDegreeOfSaturation()
    saturatedConductivityMetrePerHour = self.saturatedConductivityMetrePerDay / 24.0
    #unsaturatedConductivity = saturatedConductivityMetrePerHour * \
    #                          self.degreeOfSaturation
    b = pcr.scalar(7.0)
    self.unsaturatedConductivity = (self.degreeOfSaturation ** (2.0 * b + 3.0)) * saturatedConductivityMetrePerHour
    return self.unsaturatedConductivity, self.degreeOfSaturation

  def potentialPercolation(self):
    # calculates the percolation following van Beek 2009; potential percolation is calculated
    # which is the flux that can leave the layer to its lower layer not taking into account
    # the state of the lower layer; it accounts for moisture conditions of the current layer, i.e.
    # it will not get undersaturated if potentialPercolation is percolating (not below field capacity)
    # as percolation is limited by field capacity it will basically never get lower than field
    # capacity due to percolation
    # residual value in van Beek 2009 is moisture content at field capacity (assumed)
    potentialPercolationNotLimitedByMoistureContentFlux,tmp = self.getUnsaturatedConductivity()
    potentialPercolationNotLimitedByMoistureContentAmount = self.fluxToAmount( \
                                             potentialPercolationNotLimitedByMoistureContentFlux)
    self.calculateMaximumAbstractionThickUpToFieldCapacity()
    self.potentialPercolationAmount = pcr.min(potentialPercolationNotLimitedByMoistureContentAmount,
                                 self.maximumAbstractionThickUpToFieldCapacity)
    return self.potentialPercolationAmount

  def updateGroundWaterDepthBelowTopOfLayer(self):
    # calculates thickness of soil moisture assuming porespace above 'groundwater' (saturated zone)
    # is completely dry and then the distance (m) between top of layer and the 'groundwater' surface
    # typically used for the lower groundwater layer
    groundWaterLevelAboveBottomOfLayer = self.soilMoistureThick / self.soilPorosityFraction
    self.groundWaterDepthBelowSoil = pcr.max(0.0, self.regolithThickness - groundWaterLevelAboveBottomOfLayer)

  def updatePotentialCapillaryRiseFlux(self,unsaturatedConductivitySoilLayerMetrePerHour, saturationDegreeSoilLayer):
    # equation 9 in van Beek 2009, modified i.e. 0.5 * f is changed into a linear function returning 0 up to
    # 1 as linear function of groundwater depth between 0 and maxDepth, typically 5 m
    # depth (m) of groundwater below which no cap. rise occurs anymore
    # Takes into account water that can be substracted from groundwater (up to field capacity)
    # it is unclear why there is 0.5 in eq 9, tuned?
    # unlike van Beek 2009, the unsat conductivity of the groundwater layer is used for the conductivity
    maxDepth = pcr.scalar(5.0)
    self.updateGroundWaterDepthBelowTopOfLayer()
    # unsaturated conductivity of current layer (groundwater)
    tmp = self.getUnsaturatedConductivity()
    oneMinusSaturationDegreeSoil = pcr.max(pcr.min(1.0 - saturationDegreeSoilLayer,1.0),0.0)
    groundwaterDepthProportion = pcr.max(maxDepth - self.groundWaterDepthBelowSoil,0.0)/maxDepth
    ## formulation by me, most of the cap rise is in the groundwater so not strange to use the unsat of this layer
    #self.potentialCapillaryRiseNotCheckedFlux = pcr.scalar(self.unsaturatedConductivity) * oneMinusSaturationDegreeSoil * groundwaterDepthProportion
    # formulation following van Beek eq 9 (causes cap rise to go up if soil moisture goes up due to increase in unsat
    # but using unsat conductivity of groundwater layer
    # conductivity of the soil layer
    conductivityMetrePerHour = pcr.sqrt(self.unsaturatedConductivity * unsaturatedConductivitySoilLayerMetrePerHour)
    self.potentialCapillaryRiseNotCheckedFlux = conductivityMetrePerHour * oneMinusSaturationDegreeSoil * groundwaterDepthProportion

  def getPotentialCapillaryRiseAmount(self,unsaturatedConductivitySoilLayerMetrePerHour, saturationDegreeSoilLayer):
    # Takes into account water that can be substracted from groundwater (up to field capacity)
    self.updatePotentialCapillaryRiseFlux(unsaturatedConductivitySoilLayerMetrePerHour, saturationDegreeSoilLayer)
    self.potentialCapillaryRiseAmountNotChecked = self.fluxToAmount(self.potentialCapillaryRiseNotCheckedFlux)
    amountAboveFieldCapacity = pcr.max(self.soilMoistureThick - self.fieldCapacityThick,0.0)
    self.potentialCapillaryRiseAmount = pcr.min(amountAboveFieldCapacity,self.potentialCapillaryRiseAmountNotChecked)
    return self.potentialCapillaryRiseAmount



#  def report(self, sample, timestep, nrOfTimeSteps):
#    self.actualAdditionFlux=self.amountToFlux(self.actualAdditionFluxAmount)
#    self.lateralFlowFluxCubicMetresPerHour=self.amountToFlux(self.lateralFlowFluxAmount)*pcr.cellarea()

  def budgetCheck(self, sample, timestep):
    # NOTE this is only valid if addition,subtraction and lateral flow are invoked EACH TIME STEP
    # NOTE this does not include amountOfMoistureThickNetAdded in case of regolith thickness change
    self.actualAdditionCum = self.actualAdditionCum + self.actualAdditionFlux * self.timeStepDuration
    self.actualAbstractionCum = self.actualAbstractionCum + self.actualAbstractionFlux * self.timeStepDuration
    self.upwardSeepageCum = self.upwardSeepageCum + self.upwardSeepageFlux * self.timeStepDuration
    self.lateralFlowFluxAmountCum = self.lateralFlowFluxAmountCum + self.lateralFlowFluxAmount
    self.increaseInSubsurfaceStorage = self.soilMoistureThick - self.initialSoilMoistureThick
    budget = pcr.catchmenttotal(self.actualAdditionCum - self.increaseInSubsurfaceStorage - self.actualAbstractionCum
           - self.upwardSeepageCum, self.ldd) \
           - self.lateralFlowFluxAmountCum
    pcr.report(budget, pcrfw.generateNameST('B-sub', sample, timestep))
    pcr.report(budget / self.lateralFlowFluxAmountCum, pcrfw.generateNameST('BR-sub', sample, timestep))
    return self.increaseInSubsurfaceStorage, self.lateralFlowFluxAmountCum, self.actualAbstractionCum

#####
# adjusting regolith thickness
#####

  def possiblyAdjustSoilMoistureThickWhenRegolithThicknessChanges(self):
    amountOfMoistureThickAdded = pcr.max(0, self.wiltingPointThick - self.soilMoistureThick)
    amountOfMoistureThickRemoved = pcr.max(0, self.soilMoistureThick - self.soilPorosityThick)
    amountOfMoistureThickNetAdded = amountOfMoistureThickAdded - amountOfMoistureThickRemoved
    self.soilMoistureThick = pcr.min(pcr.max(self.soilMoistureThick, self.wiltingPointThick), self.soilPorosityThick)
    return amountOfMoistureThickNetAdded

  def updateRegolithThickness(self, newRegolithThickness):
     # Updates the state of the regolith 'aiming at' a situation where the amount of water in the soil
     # remains the same, it is only removed when it is above the porosity and only added when it is below the wilting
     # point (both physically impossible soil moisture values) but otherwise it remains the same amount of
     # water in the subsoil.
     self.convertFractionsToThickness(newRegolithThickness, self.soilPorosityFraction, self.wiltingPointFraction,
                                  self.limitingPointFraction, self.fieldCapacityFraction)
     amountOfMoistureThickNetAdded = self.possiblyAdjustSoilMoistureThickWhenRegolithThicknessChanges()
     return amountOfMoistureThickNetAdded

  def updateBedrock(self, newLdd, newDemOfBedrock):
    self.setBedrock(newLdd, newDemOfBedrockTopography)

###
# calculations for additional reports

  def calculateSoilMoistureFraction(self):
    self.soilMoistureFraction = self.soilMoistureThick / self.regolithThickness



# ldd='ldd.map'
# demOfBedrockTopography='mdtpaz4.map'
# regolithThickness=pcr.scalar(10.0)
# initialSoilMoistureFraction=pcr.scalar(0.5)
# soilPorosityFraction=pcr.scalar(0.6)
# wiltingPointFraction=pcr.scalar(0.2)
# fieldCapacityFraction=pcr.scalar(0.3)
# saturatedConductivityMetrePerDay=pcr.scalar(2.0)
# timeStepDuration=1
#
# d_surfaceWaterOneLayer=SubsurfaceWaterOneLayer(
#              ldd,
#              demOfBedrockTopography,
#              regolithThickness,
#              initialSoilMoistureFraction,
#              soilPorosityFraction,
#              wiltingPointFraction,
#              fieldCapacityFraction,
#              saturatedConductivityMetrePerDay,
#              timeStepDuration)

# test on addition of water
# maximumAdditionFlux=d_surfaceWaterOneLayer.getMaximumAdditionFlux()
# report(maximumAdditionFlux,'maf.map')
# actualAdditionFlux=d_surfaceWaterOneLayer.addWater(0.99)
# report(actualAdditionFlux,'aaf.map')
#
# test on abstraction of water
# maximumAbstractionFlux=d_surfaceWaterOneLayer.getMaximumAbstractionFlux()
# report(maximumAbstractionFlux,'mabf.map')
# actualAbstractionFlux=d_surfaceWaterOneLayer.abstractWater(2.99)
# report(actualAbstractionFlux,'aabf.map')


# class testModel():
#  def __init__(self):
#    pass
#
#  def premcloop(self):
#    pass
#
#  def initial(self):
#    self.ldd='ldd.map'
#    demOfBedrockTopography='mdtpaz4.map'
#    regolithThickness=pcr.scalar(2.0)
#    initialSoilMoistureFraction=pcr.scalar(0.5)
#    soilPorosityFraction=pcr.scalar(0.6)
#    wiltingPointFraction=pcr.scalar(0.2)
#    fieldCapacityFraction=pcr.scalar(0.3)
#    saturatedConductivityMetrePerDay=pcr.scalar(100.0)
#    timeStepDuration=1
#
#    self.d_subsurfaceWaterOneLayer=SubsurfaceWaterOneLayer(
#                  self.ldd,
#                  demOfBedrockTopography,
#                  regolithThickness,
#                  initialSoilMoistureFraction,
#                  soilPorosityFraction,
#                  wiltingPointFraction,
#                  fieldCapacityFraction,
#                  saturatedConductivityMetrePerDay,
#                  timeStepDuration)
#
#  def dynamic(self):
#    maximumAdditionFlux=self.d_subsurfaceWaterOneLayer.getMaximumAdditionFlux()
#    self.report(maximumAdditionFlux,'ma')
#    actualAdditionFlux=self.d_subsurfaceWaterOneLayer.addWater(ifthenelse(pcrlt(mapuniform(),0.1),pcr.scalar(0.1),0.0))
#    actualAbstractionFlux=self.d_subsurfaceWaterOneLayer.abstractWater(0.011)
#    upwardSeepageFlux=self.d_subsurfaceWaterOneLayer.lateralFlow()
#    self.d_subsurfaceWaterOneLayer.report(self.currentSampleNumber(), self.currentTimeStep())
#    self.d_subsurfaceWaterOneLayer.budgetCheck(self.currentSampleNumber() , self.currentTimeStep())
#    self.report(pcr.accuflux(self.ldd,upwardSeepageFlux),'q')
#
#  def postmcloop(self):
#    pass
#
#myModel = testModel()
#dynamicModel = DynamicFramework(myModel, 24*30)
#mcModel = MonteCarloFramework(dynamicModel, 1)
# mcModel.run()
