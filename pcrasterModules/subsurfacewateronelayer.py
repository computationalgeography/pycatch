import math
from pcraster import *
from pcraster.framework import *
import generalfunctions, component

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


setclone('clone.map')

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
              timeStepDuration,
              timeStepsToReport,
              setOfVariablesToReport):

    ######
    # inits only for suspend and resume in filtering
    ######

    self.actualAbstractionFlux=scalar(0)
    self.actualAbstractionFluxAmount=scalar(0)
    self.actualAdditionFlux=scalar(0)
    self.actualAdditionFluxAmount=scalar(0)
    self.fWaterPotential=scalar(0)
    self.lateralFlowCubicMetrePerDay=scalar(0)
    self.lateralFlowCubicMetrePerTimeStep=scalar(0)
    self.lateralFlowDarcyFluxAmount=scalar(0)
    self.lateralFlowFluxAmount=scalar(0)
    self.lateralFlowFluxCubicMetresPerHour=scalar(0)
    self.maximumAbstractionThick=scalar(0)
    self.maximumAdditionThick=scalar(0)
    self.saturatedLayerThick=scalar(0)
    self.saturatedLayerThickness=scalar(0)
    self.soilMoistureMinusWiltingPointThick=scalar(0)
    self.upwardSeepageAmount=scalar(0)
    self.upwardSeepageFlux=scalar(0)
    self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour=scalar(0)
    self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour=scalar(0)
    self.totalSoilMoistureThickInUpstreamAreaCubicMetre=scalar(0)
    self.totalSoilMoistureFractionInUpstreamAreaTimesCellArea=scalar(0)
    self.totalSaturatedThickInUpstreamAreaCubicMetre=scalar(0)
    self.soilMoistureFraction=scalar(0)
    self.variablesToReport={}
    self.variablesAsNumpyToReport={}
    self.wiltingPointThick=scalar(0)

    ####
    # real inits
    ####

    self.cellArea=cellarea()

    self.setBedrock(ldd,demOfBedrockTopography)

    # parameters depending on regolith thickness
    self.convertFractionsToThickness(regolithThickness,soilPorosityFraction,wiltingPointFraction, \
                                  limitingPointFraction, fieldCapacityFraction)

    # main state variable, amount of water
    self.initialSoilMoistureFraction=initialSoilMoistureFraction
    self.initialSoilMoistureThick=self.initialSoilMoistureFraction*self.regolithThickness
    self.soilMoistureThick=self.initialSoilMoistureThick

    # values independent of regolith thickness 
    self.saturatedConductivityMetrePerDay=saturatedConductivityMetrePerDay
    self.timeStepDuration=scalar(timeStepDuration)
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport

    # DJ add check on fractions (e.g. wiltingpoint cannot be bigger than soil porosity)

    # for budgets
    self.actualAdditionCum=scalar(0)
    self.actualAbstractionCum=scalar(0)
    self.lateralFlowFluxAmountCum=scalar(0)
    self.upwardSeepageCum=scalar(0)


  def reportAsMaps(self,sample,timestep):
    # reports
    self.variablesToReport = {}
    if self.setOfVariablesToReport == 'full':
      self.variablesToReport = {
                                'Gs': self.soilMoistureThick,
                                #'Gad': self.maximumAdditionThick,
                                'Go': self.actualAbstractionFlux,  
                                #'Gi': self.actualAdditionFlux,    
                                #'Gq': self.lateralFlowFluxCubicMetresPerHour,
                                #'Gwp': self.fWaterPotential,
                                #'rts': self.regolithThickness,
                                #'Gks': self.saturatedConductivityMetrePerDay
                                'Gxt': self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour,
                                'Got': self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour
                                 }
    if self.setOfVariablesToReport == 'filtering':
      self.variablesToReport = { }
    self.reportMaps(sample,timestep)

  def updateVariablesAsNumpyToReport(self):
    self.variablesAsNumpyToReport = {
                                'Gxt': self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour,
                                'Got': self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour,
                                'Gst': self.totalSoilMoistureThickInUpstreamAreaCubicMetre,
                                'Gtt': self.totalSoilMoistureFractionInUpstreamAreaTimesCellArea,
                                'GSt': self.totalSaturatedThickInUpstreamAreaCubicMetre
                                    }

  def reportAsNumpyOneFile(self,locations,sample,timestep,endTimeStep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpyOneFilePerRealization(locations,sample,timestep,endTimeStep)

  def reportAsNumpyMultipleFiles(self,locations,sample,timestep):
    self.updateVariablesAsNumpyToReport()
    self.reportAsNumpy(locations,sample,timestep)


  def convertFractionsToThickness(self,regolithThickness,soilPorosityFraction,wiltingPointFraction, \
                                  limitingPointFraction, fieldCapacityFraction):
    # regolith thickness and parameters derived from regolith thickness,
    # need to be recalculated when regolith thickness changes
    self.regolithThickness=regolithThickness

    self.soilPorosityFraction=soilPorosityFraction
    self.soilPorosityThick=self.soilPorosityFraction*self.regolithThickness

    self.wiltingPointFraction=wiltingPointFraction
    self.wiltingPointThick=scalar(self.wiltingPointFraction*self.regolithThickness)

    self.limitingPointFraction=limitingPointFraction
    self.limitingPointThick=self.limitingPointFraction*self.regolithThickness

    self.fieldCapacityFraction=fieldCapacityFraction
    self.fieldCapacityThick=self.fieldCapacityFraction*self.regolithThickness

    self.maximumAvailableForVegetationUptakeThick=self.limitingPointThick-self.wiltingPointThick

  def setBedrock(self,ldd,demOfBedrockTopography):
    self.ldd=ldd
    self.demOfBedrockTopography=demOfBedrockTopography
    self.slopeToDownstreamNeighbourNotFlat=generalfunctions.slopeToDownstreamNeighbourNotFlat( \
                                           self.demOfBedrockTopography,self.ldd,0.001)


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

  def fluxToAmount(self,flux):
    fluxAmount=flux * self.timeStepDuration
    return fluxAmount

  def amountToFlux(self,amount):
    flux=amount/self.timeStepDuration
    return flux

  def thicknessToFraction(self,thickness):
    fraction=thickness/self.regolithThickness
    return fraction

  def calculateMaximumAbstractionThick(self):
    # amount of water that can be taken out of the storage, up to wilting point
    self.maximumAbstractionThick=self.soilMoistureThick-self.wiltingPointThick

  def getMaximumAbstractionFlux(self):
    self.calculateMaximumAbstractionThick()
    maximumAbstractionFlux=self.amountToFlux(self.maximumAbstractionThick)
    return maximumAbstractionFlux

  def calculateMaximumAdditionThick(self):
    # amount of water that can be added, up to saturation
    self.maximumAdditionThick=self.soilPorosityThick-self.soilMoistureThick

  def getMaximumAdditionFlux(self):
    self.calculateMaximumAdditionThick()
    maximumAdditionFlux=self.amountToFlux(self.maximumAdditionThick)
    return maximumAdditionFlux

  def abstractWater(self,potentialAbstractionFlux):
    # DJ add check on potentialAbstractionFlux >= 0
    self.calculateMaximumAbstractionThick()
    potentialAbstractionFluxAmount=self.fluxToAmount(potentialAbstractionFlux)
    self.actualAbstractionFluxAmount=min(self.maximumAbstractionThick,potentialAbstractionFluxAmount)
    self.soilMoistureThick=self.soilMoistureThick-self.actualAbstractionFluxAmount

    # conversions
    self.actualAbstractionFlux=self.amountToFlux(self.actualAbstractionFluxAmount)

    # for a report
    self.totalActualAbstractionInUpstreamArea()

    return self.actualAbstractionFlux

  def addWater(self,potentialAdditionFlux):
    # DJ add check on potentialAdditionFlux >= 0
    self.calculateMaximumAdditionThick()
    potentialAdditionFluxAmount=self.fluxToAmount(potentialAdditionFlux)
    self.actualAdditionFluxAmount=min(self.maximumAdditionThick,potentialAdditionFluxAmount)
    self.soilMoistureThick=self.soilMoistureThick+self.actualAdditionFluxAmount

    # conversions
    self.actualAdditionFlux=self.amountToFlux(self.actualAdditionFluxAmount)

    return self.actualAdditionFlux

  def calculateThicknessOfSaturatedLayer(self):
    # thickness of saturated layer in m water depth (no pores)
    self.saturatedLayerThick=max(self.soilMoistureThick-self.fieldCapacityFraction*self.regolithThickness,scalar(0.0))
    # thickness of saturated layer in m (including pores)
    self.saturatedLayerThickness=self.saturatedLayerThick/self.soilPorosityFraction

  def calculateLateralFlowDarcy(self):
    self.calculateThicknessOfSaturatedLayer()
    self.lateralFlowCubicMetrePerDay=self.slopeToDownstreamNeighbourNotFlat*self.saturatedConductivityMetrePerDay \
                                * celllength() * self.saturatedLayerThickness;
    self.lateralFlowCubicMetrePerTimeStep=(self.lateralFlowCubicMetrePerDay/24.0)*self.timeStepDuration
    self.lateralFlowDarcyFluxAmount=self.lateralFlowCubicMetrePerTimeStep/cellarea()

  def calculateSoilMoistureMinusWiltingPointThick(self):
    self.soilMoistureMinusWiltingPointThick=max(self.soilMoistureThick-self.wiltingPointThick,scalar(0.0))

  def calculateLateralFlow(self):
    self.calculateLateralFlowDarcy()
    # outgoing lateral flow may not be greater than soil moisture minus wilting point
    self.calculateSoilMoistureMinusWiltingPointThick()
    self.lateralFlowFluxAmount=min(self.lateralFlowDarcyFluxAmount,self.soilMoistureMinusWiltingPointThick)
    self.lateralFlowFluxCubicMetresPerHour=self.amountToFlux(self.lateralFlowFluxAmount)*cellarea()
    
  def lateralFlow(self):
    self.calculateLateralFlow()
    self.soilMoistureThick=self.soilMoistureThick-self.lateralFlowFluxAmount+upstream(self.ldd,self.lateralFlowFluxAmount)
    self.upwardSeepageAmount=max(scalar(0.0),self.soilMoistureThick-self.soilPorosityThick)
    self.soilMoistureThick=min(self.soilMoistureThick,self.soilPorosityThick)

    self.upwardSeepageFlux=self.amountToFlux(self.upwardSeepageAmount)
    self.totalUpwardSeepageInUpstreamArea()
    self.totalSoilMoistureThickInUpstreamArea()
    self.totalSoilMoistureFractionInUpstreamArea()
    self.totalSaturatedThickInUpstreamArea()
    return self.upwardSeepageFlux

  def totalUpwardSeepageInUpstreamArea(self):
    self.totalUpwardSeepageInUpstreamAreaCubicMetrePerHour=accuflux(self.ldd, self.upwardSeepageFlux)*self.cellArea

  def totalActualAbstractionInUpstreamArea(self):
    self.totalActualAbstractionInUpstreamAreaCubicMetrePerHour=accuflux(self.ldd, self.actualAbstractionFlux)*self.cellArea

  def totalSoilMoistureThickInUpstreamArea(self):
    self.totalSoilMoistureThickInUpstreamAreaCubicMetre=accuflux(self.ldd, self.soilMoistureThick)*self.cellArea

  def totalSoilMoistureFractionInUpstreamArea(self):
    self.calculateSoilMoistureFraction()
    self.totalSoilMoistureFractionInUpstreamAreaTimesCellArea=accuflux(self.ldd, self.soilMoistureFraction)*self.cellArea

  def totalSaturatedThickInUpstreamArea(self):
    self.totalSaturatedThickInUpstreamAreaCubicMetre=accuflux(self.ldd, self.saturatedLayerThick)*self.cellArea
  
  def getFWaterPotential(self):
    # calculates the reduction factor (0-1) for the stomatal conductance, i.e. stomatal
    # conductance in Penman will be the maximum stomatal conductance multiplied by
    # this reduction factor (see penman script)
    # the reduction factor is a linear function of soil moisture content, zero at
    # or below wilting point and one at or above field capacity, and a linear function inbetween
    fWaterPotentialTmp=(self.soilMoistureThick-self.wiltingPointThick) \
                           / self.maximumAvailableForVegetationUptakeThick
    fWaterPotentialTmpWilt=ifthenelse(pcrlt(self.soilMoistureThick,self.wiltingPointThick), \
                           scalar(0),fWaterPotentialTmp)
    self.fWaterPotential=ifthenelse(pcrgt(self.soilMoistureThick,self.limitingPointThick), \
                           scalar(1),fWaterPotentialTmpWilt)
    return self.fWaterPotential

#  def report(self, sample, timestep, nrOfTimeSteps):
#    self.actualAdditionFlux=self.amountToFlux(self.actualAdditionFluxAmount)
#    self.lateralFlowFluxCubicMetresPerHour=self.amountToFlux(self.lateralFlowFluxAmount)*cellarea()

  def budgetCheck(self, sample, timestep):
    # NOTE this is only valid if addition,subtraction and lateral flow are invoked EACH TIME STEP
    # NOTE this does not include amountOfMoistureThickNetAdded in case of regolith thickness change
    self.actualAdditionCum=self.actualAdditionCum+self.actualAdditionFlux*self.timeStepDuration
    self.actualAbstractionCum=self.actualAbstractionCum+self.actualAbstractionFlux*self.timeStepDuration
    self.upwardSeepageCum=self.upwardSeepageCum+self.upwardSeepageFlux*self.timeStepDuration
    self.lateralFlowFluxAmountCum=self.lateralFlowFluxAmountCum+self.lateralFlowFluxAmount
    self.increaseInSubsurfaceStorage=self.soilMoistureThick-self.initialSoilMoistureThick
    budget=catchmenttotal(self.actualAdditionCum-self.increaseInSubsurfaceStorage-self.actualAbstractionCum \
           -self.upwardSeepageCum,self.ldd) \
           -self.lateralFlowFluxAmountCum
    report(budget,generateNameST('B-sub', sample, timestep))
    report(budget/self.lateralFlowFluxAmountCum,generateNameST('BR-sub', sample, timestep))
    return self.increaseInSubsurfaceStorage, self.lateralFlowFluxAmountCum, self.actualAbstractionCum

#####
# adjusting regolith thickness    
#####

  def possiblyAdjustSoilMoistureThickWhenRegolithThicknessChanges(self):
    amountOfMoistureThickAdded=max(0,self.wiltingPointThick-self.soilMoistureThick)
    amountOfMoistureThickRemoved=max(0,self.soilMoistureThick-self.soilPorosityThick)
    amountOfMoistureThickNetAdded=amountOfMoistureThickAdded-amountOfMoistureThickRemoved
    self.soilMoistureThick=min(max(self.soilMoistureThick,self.wiltingPointThick),self.soilPorosityThick)
    return amountOfMoistureThickNetAdded
  
  def updateRegolithThickness(self,newRegolithThickness):
     self.convertFractionsToThickness(newRegolithThickness,self.soilPorosityFraction,self.wiltingPointFraction, \
                                  self.limitingPointFraction, self.fieldCapacityFraction)
     amountOfMoistureThickNetAdded=self.possiblyAdjustSoilMoistureThickWhenRegolithThicknessChanges()
     return amountOfMoistureThickNetAdded

  def updateBedrock(self,newLdd,newDemOfBedrock):
    self.setBedrock(newLdd,newDemOfBedrockTopography)

###
# calculations for additional reports

  def calculateSoilMoistureFraction(self):
    self.soilMoistureFraction=self.soilMoistureThick/self.regolithThickness
    

    
#ldd='ldd.map'
#demOfBedrockTopography='mdtpaz4.map'
#regolithThickness=scalar(10.0)
#initialSoilMoistureFraction=scalar(0.5)
#soilPorosityFraction=scalar(0.6)
#wiltingPointFraction=scalar(0.2)
#fieldCapacityFraction=scalar(0.3)
#saturatedConductivityMetrePerDay=scalar(2.0)
#timeStepDuration=1
#
#d_surfaceWaterOneLayer=SubsurfaceWaterOneLayer(
#              ldd,
#              demOfBedrockTopography,
#              regolithThickness,
#              initialSoilMoistureFraction,
#              soilPorosityFraction,
#              wiltingPointFraction,
#              fieldCapacityFraction, 
#              saturatedConductivityMetrePerDay,
#              timeStepDuration)

## test on addition of water
#maximumAdditionFlux=d_surfaceWaterOneLayer.getMaximumAdditionFlux()
#report(maximumAdditionFlux,'maf.map')
#actualAdditionFlux=d_surfaceWaterOneLayer.addWater(0.99)
#report(actualAdditionFlux,'aaf.map')
#
## test on abstraction of water
#maximumAbstractionFlux=d_surfaceWaterOneLayer.getMaximumAbstractionFlux()
#report(maximumAbstractionFlux,'mabf.map')
#actualAbstractionFlux=d_surfaceWaterOneLayer.abstractWater(2.99)
#report(actualAbstractionFlux,'aabf.map')


#class testModel():
#  def __init__(self):
#    pass
#
#  def premcloop(self):
#    pass
#
#  def initial(self):
#    self.ldd='ldd.map'
#    demOfBedrockTopography='mdtpaz4.map'
#    regolithThickness=scalar(2.0)
#    initialSoilMoistureFraction=scalar(0.5)
#    soilPorosityFraction=scalar(0.6)
#    wiltingPointFraction=scalar(0.2)
#    fieldCapacityFraction=scalar(0.3)
#    saturatedConductivityMetrePerDay=scalar(100.0)
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
#    actualAdditionFlux=self.d_subsurfaceWaterOneLayer.addWater(ifthenelse(pcrlt(mapuniform(),0.1),scalar(0.1),0.0))
#    actualAbstractionFlux=self.d_subsurfaceWaterOneLayer.abstractWater(0.011)
#    upwardSeepageFlux=self.d_subsurfaceWaterOneLayer.lateralFlow()
#    self.d_subsurfaceWaterOneLayer.report(self.currentSampleNumber(), self.currentTimeStep())
#    self.d_subsurfaceWaterOneLayer.budgetCheck(self.currentSampleNumber() , self.currentTimeStep())
#    self.report(accuflux(self.ldd,upwardSeepageFlux),'q')
#
#  def postmcloop(self):
#    pass
#
#myModel = testModel()
#dynamicModel = DynamicFramework(myModel, 24*30)
#mcModel = MonteCarloFramework(dynamicModel, 1)
#mcModel.run()
