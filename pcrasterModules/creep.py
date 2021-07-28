from pcraster import *
import sys, generalfunctions
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

class Creep(component.Component):
  def __init__(self,dem,timeStepDuration,diffusion, timeStepsToReport,setOfVariablesToReport):
    import generalfunctions
    self.dem=dem
    self.timeStepDuration=timeStepDuration
    self.diffusion=diffusion
    self.timeStepsToReport=timeStepsToReport
    self.setOfVariablesToReport=setOfVariablesToReport
    self.upperNBMissing,self.rightNBMissing,self.lowerNBMissing,self.leftNBMissing= \
                   generalfunctions.neighbourIsMissingValueOrEdgeAndCellItselfIsDefined(self.dem)
    self.soilThickness=scalar(0)
    self.outflow=scalar(0)
    self.flowOverBoundaries=scalar(0)
    self.correctedFactor=scalar(0)
    self.velocityMetrePerYear=scalar(0)

    self.output_mapping = {
                           'Ds': self.soilThickness,
                           'Dou': self.outflow,
                           'Dbo': self.flowOverBoundaries,
                           'Dco': self.correctedFactor,
                           'Dve': self.velocityMetrePerYear
                             }

  def reportAsMaps(self, sample, timestep):
    self.variablesToReport = self.rasters_to_report(self.setOfVariablesToReport)
    self.reportMaps(sample, timestep)


  def gradient(self):
    # slopeX, transport to the left is positive
    # slopeY, transport to the bottom is positive
    slopeX=gradx(self.dem)
    slopeY=grady(self.dem)
    return slopeX, slopeY

  def amount(self,soilThickness):
    slopeX,slopeY=self.gradient() 
    amountX=self.diffusion*slopeX*soilThickness*self.timeStepDuration
    amountY=self.diffusion*slopeY*soilThickness*self.timeStepDuration
    return amountX,amountY, slopeX, slopeY

  def amountNotGreaterThanFractionOfSoilThickness(self,soilThickness,amountX,amountY):
    totAmount=abs(amountX)+abs(amountY)
    #totAmount=amountX+amountY
    maxAllowed=soilThickness*0.5
    correctionFactor=maxAllowed/totAmount
    amountXCorrectedInCase=correctionFactor*amountX
    amountYCorrectedInCase=correctionFactor*amountY
    doCorrection=totAmount > maxAllowed
    amountXCorrected=ifthenelse(doCorrection,amountXCorrectedInCase,amountX)
    amountYCorrected=ifthenelse(doCorrection,amountYCorrectedInCase,amountY)
    correctedFactor=ifthenelse(doCorrection,correctionFactor,scalar(1))
    return amountXCorrected, amountYCorrected, correctedFactor

  def westingInShift(self,amountX):
    westing=boolean(amountX > 0.0)
    return westing

  def northingInShift(self,amountY):
    northing=boolean(amountY < 0.0)
    return northing

  def shiftWesting(self,amount,westing):
    # move to left
    removedWestingTrue=ifthen(westing,amount)
    # move to right
    removedWestingFalse=ifthen(pcrnot(westing),-amount)
    transportedWestingTrue=shift0(removedWestingTrue,scalar(0),1.0)
    transportedWestingFalse=shift0(removedWestingFalse,scalar(0),-1.0)
    return transportedWestingTrue+transportedWestingFalse

  def shiftNorthing(self,amount,northing):
    removedNorthingTrue=ifthen(northing,-amount)
    removedNorthingFalse=ifthen(pcrnot(northing),amount)
    transportedNorthingTrue=shift0(removedNorthingTrue,1.0,scalar(0))
    transportedNorthingFalse=shift0(removedNorthingFalse,-1.0,scalar(0))
    return transportedNorthingTrue+transportedNorthingFalse

  def transportX(self,amountX):
    westing=self.westingInShift(amountX)
    inflow=self.shiftWesting(amountX,westing)
    return inflow

  def transportY(self,amountY):
    northing=self.northingInShift(amountY)
    inflow=self.shiftNorthing(amountY,northing)
    return inflow

  def advectionOutflow(self,amountX,amountY,upperNBMissing,rightNBMissing,lowerNBMissing,leftNBMissing):
    amountToUpperNB=ifthenelse(pcrand(upperNBMissing,amountY < 0.0),0-amountY,scalar(0))
    amountToLowerNB=ifthenelse(pcrand(lowerNBMissing,amountY > 0.0),amountY,scalar(0))
    amountToRightNB=ifthenelse(pcrand(rightNBMissing,amountX < 0.0),0-amountX,scalar(0))
    amountToLeftNB=ifthenelse(pcrand(leftNBMissing,amountX > 0.0),amountX,scalar(0))
    return amountToUpperNB + amountToLowerNB + amountToRightNB + amountToLeftNB

  def advection(self,material,amountX,amountY):
    inflowX=self.transportX(amountX)
    inflowY=self.transportY(amountY)
    #material=max(material-abs(amountX)-abs(amountY)+inflowX+inflowY,0.0)
    material=max(material-abs(amountX)-abs(amountY)+inflowX+inflowY,0.0)
    return material, inflowX, inflowY

  def setFlowOverBoundariesToZero(self,amountX,amountY):
    amountY=ifthenelse(self.upperNBMissing,max(amountY,scalar(0)),amountY)
    amountY=ifthenelse(self.lowerNBMissing,min(amountY,scalar(0)),amountY)
    amountX=ifthenelse(self.rightNBMissing,max(amountX,scalar(0)),amountX)
    amountX=ifthenelse(self.leftNBMissing, min(amountX,scalar(0)),amountX)
    return amountX, amountY

  def diffuse(self,soilThickness,bedrockDem,closedBoundaries):
    self.dem=bedrockDem
    amountXUncor,amountYUncor,slopeX,slopeY=self.amount(soilThickness)
    amountX,amountY,self.correctedFactor=self.amountNotGreaterThanFractionOfSoilThickness(soilThickness,amountXUncor,amountYUncor)
    if closedBoundaries:
      amountX,amountY=self.setFlowOverBoundariesToZero(amountX,amountY)
    self.soilThickness,inflowX,inflowY=self.advection(soilThickness,amountX,amountY)
    self.flowOverBoundaries=self.advectionOutflow(amountX,amountY, \
                       self.upperNBMissing,self.rightNBMissing,self.lowerNBMissing,self.leftNBMissing)
    self.outflow=abs(amountX)+abs(amountY)
    self.velocityMetrePerYear=((self.outflow*celllength())/soilThickness)/self.timeStepDuration
    return self.soilThickness, self.outflow, self.flowOverBoundaries, self.correctedFactor, amountX, amountY, inflowX,inflowY
   

### test
#dem='mdtpaz4.map'
#soilThick=ifthen(defined(dem),uniform(1))
#timeStepDuration=1.0
#report(soilThick,'thickb.map')
#
#d_creep = Creep(dem,timeStepDuration,0.1,'a','b')
#soilThick,outflow,flowOverBoundaries,correctedFactor, amountX,amountY,inflowX,inflowY=d_creep.diffuse(soilThick,dem)
#report(amountX,'amountX.map')
#report(amountY,'amountY.map')
#report(inflowX,'inflowX.map')
#report(inflowY,'inflowY.map')
#report(soilThick,'thicka.map')
#
#report(outflow,'out.map')
#report(correctedFactor,'corr.map')
#
#totAfter=maptotal(soilThick+flowOverBoundaries)
#report(totAfter,'totaft.map')

#slopeX,slopeY=d_creep.gradient()
#report(slopeX,'sx.map')
#report(slopeY,'sy.map')

#westing=d_creep.westingInShift(slopeX)
#northing=d_creep.northingInShift(slopeY)
#uniForm=uniform(ifthen(defined(dem),boolean(1)))
#report(uniForm,'uni.map')
#report(northing,'northing.map')
#piet=ifthen(northing,uniForm)
#removedNorthingTrue, removedNorthingFalse,transportedNorthingTrue, transportedNorthingFalse \
#        = d_creep.shiftNorthing(uniform,northing)
#report(remNorthingTrue,'rwt.map')
#report(remNorthingFalse,'rwf.map')
#report(transpNorthingTrue,'twt.map')
#report(transpNorthingFalse,'twf.map')

