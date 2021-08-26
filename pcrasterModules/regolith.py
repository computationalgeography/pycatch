from pcraster import * 
import sys
from pcraster.framework import *
import component

# notes
# time step duration in years
# levels in m
# vertical fluxes in m/year, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a layer of soil (i.e. air + rock material over a whole cell)
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

class RegolithDemAndBedrock(component.Component):
  def __init__(self, dem, regolithThickness, timeStepDuration, timeStepsToReport,setOfVariablesToReport):

    # real inits
    self.dem=dem
    self.regolithThickness=scalar(regolithThickness)
    self.demOfBedrock=self.dem-self.regolithThickness
    self.timeStepDuration=timeStepDuration
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport
    #self.minimumAllowedRegolithThickness=0.000001
    self.minimumAllowedRegolithThickness=0.001
    self.surfaceLdd=lddcreate(self.dem,1e31,1e31,1e31,1e31)
    self.bedrockLdd=lddcreate(self.demOfBedrock,1e31,1e31,1e31,1e31)

    # for reports
    self.initialDem=self.dem

    # budget check
    # needs to be written

  def reportAsMaps(self, sample, timestep):
    self.output_mapping = {
                          'Ast': self.regolithThickness,
                          'Adb': self.demOfBedrock,
                          'Ade': self.dem,
                          'Abl': self.bedrockLdd,
                          'Asl': self.surfaceLdd,
                           }
    self.variablesToReport = self.rasters_to_report(self.setOfVariablesToReport)
    self.reportMaps(sample, timestep)

  def amountToFlux(self,amount):
    flux=amount/self.timeStepDuration
    return flux

  def fluxToAmount(self,flux):
    amount=flux * self.timeStepDuration
    return amount

  def updateWithBedrockWeathering(self,bedrockWeatheringFlux):
    self.regolithThickness=self.regolithThickness+self.fluxToAmount(bedrockWeatheringFlux)
    # may 2017, added max statement for minimum allowed reg thickness
    self.regolithThickness=max(self.regolithThickness,self.minimumAllowedRegolithThickness)
    self.demOfBedrock=self.dem-self.regolithThickness
    self.bedrockLdd=lddcreate(self.demOfBedrock,1e31,1e31,1e31,1e31)

  def updateWithDeposition(self,potentialDepositionFlux):
    '''depositionFlux can be positive or negative'''
    potentialDepositionAmount=self.fluxToAmount(potentialDepositionFlux)
    actualDepositionAmount=max(self.minimumAllowedRegolithThickness-self.regolithThickness,potentialDepositionAmount)
    self.regolithThickness=self.regolithThickness+actualDepositionAmount
    self.dem=self.demOfBedrock+self.regolithThickness
    self.surfaceLdd=lddcreate(self.dem,1e31,1e31,1e31,1e31)
    return self.amountToFlux(actualDepositionAmount)
  
  def setNewRegolith(self,newRegolithThickness):
    self.regolithThickness=newRegolithThickness
    self.dem=self.demOfBedrock+self.regolithThickness
    self.surfaceLdd=lddcreate(self.dem,1e31,1e31,1e31,1e31)

  def setBaselevel(self,baseLevel):
    '''
    surface dem is adjusted to baseLevel only where baseLevel is not mv
    all other properties are adjusted, ie bedrock removal is done when required
    note that bedrockAddedAmount is always negative (i.e., always removal)
    '''
    # calculation for area with non mv's on baseLevel
    demOfBedrockBaseLevelArea=min(self.demOfBedrock,baseLevel) 
    regolithThicknessBaseLevelArea=max(baseLevel-demOfBedrockBaseLevelArea,self.minimumAllowedRegolithThickness)

    # return values
    bedrockAddedAmount=cover(demOfBedrockBaseLevelArea-self.demOfBedrock,0)
    regolithAddedAmount=cover(regolithThicknessBaseLevelArea-self.regolithThickness,0)

    # fill in mv's on baselevel with original values
    self.demOfBedrock=cover(demOfBedrockBaseLevelArea,self.demOfBedrock)
    self.regolithThickness=cover(regolithThicknessBaseLevelArea,self.regolithThickness)

    # update other state variables
    self.dem=self.demOfBedrock+self.regolithThickness
    self.surfaceLdd=lddcreate(self.dem,1e31,1e31,1e31,1e31)
    self.bedrockLdd=lddcreate(self.demOfBedrock,1e31,1e31,1e31,1e31)
    return regolithAddedAmount, bedrockAddedAmount
    
  def getRegolithProperties(self):
    return self.regolithThickness,self.demOfBedrock,self.dem, self.bedrockLdd,self.surfaceLdd

#  def budgetCheck(self, sample, timestep):
   # needs to take into account actualDepositionFlux
