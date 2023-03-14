from pcraster import *
#from genericfunctions import *
import sys
import generalfunctions

# notes
# time step duration in h
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a waterslice over a whole cell)
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types
# the module has to update self.store by itself!
# also, the module has to update cumulative discharge by itself 

def wettedPerimeter(waterHeight,channelWidth):
  wettedPerimeter=2*waterHeight+channelWidth
  return wettedPerimeter

class RunoffKinematic:
  def __init__(self, ldd, dem, channelWidth, n, beta, initialDischargeCubicMetrePerSecond, 
               initialWaterHeight, timeStepDuration):
    self.ldd=ldd
    self.dem=dem
    self.channelWidth=channelWidth
    self.n=n
    self.beta=beta
    self.Q=initialDischargeCubicMetrePerSecond
    self.waterHeightInChannel=initialWaterHeight
    self.timeStepDuration=timeStepDuration
    self.timeStepDurationSeconds=self.timeStepDuration*scalar(3600)
    self.slopeToDownstreamNeighbour=generalfunctions.slopeToDownstreamNeighbourNotFlat(self.dem,self.ldd,0.0001)
    self.distanceToDownstreamCell=generalfunctions.distancetodownstreamcell(self.ldd)
    self.alpTerm=(self.n/(sqrt(self.slopeToDownstreamNeighbour)))**self.beta;
    self.alpPow=(scalar(2)/scalar(3))*self.beta;
    self.wettedPerimeter=wettedPerimeter(self.waterHeightInChannel,self.channelWidth)
    self.alpha=self.alpTerm*(self.wettedPerimeter**self.alpPow)
    self.cellArea=cellarea()
    self.channelArea=self.channelWidth*self.distanceToDownstreamCell
    self.store=self.waterHeightInChannel*self.channelArea
    self.cumulativeDischargeCubicMetres=scalar(0.0)
    print('note that there may be an error in runonflux in kinematic, can be ignored for now')

  def runonFlux(self):
    runonFlux=self.Q/self.cellArea
    return runonFlux

  def perHourToPerSecond(self,aScalar):
    perSecond=aScalar/scalar(3600)
    return perSecond

  def perSecondToPerHour(self,aScalar):
    perHour=aScalar*scalar(3600)
    return perHour

  def perHourToPerTimeStep(self,aScalar):
    perTimeStep=aScalar*self.timeStepDuration
    return perTimeStep

  def perSecondToPerTimeStep(self,aScalar):
    perHour=self.perSecondToPerHour(aScalar)
    perTimeStep=self.perHourToPerTimeStep(perHour)
    return perTimeStep

  def calculateActualAbstractionFlux(self,rainFlux,potentialAbstractionFlux):
    runonFlux=self.runonFlux()
    #print('THIS IS WRONG UNIT, runoFlux....')
    potentialAvailableForAbstractionFlux=runonFlux + rainFlux
    actualAbstractionFlux=min(potentialAbstractionFlux,potentialAvailableForAbstractionFlux)
    return actualAbstractionFlux

  def sideFlow(self,rainFlux,potentialAbstractionFlux):
    self.actualAbstractionFlux=self.calculateActualAbstractionFlux(rainFlux,potentialAbstractionFlux)
    changeStreamFlow=rainFlux-self.actualAbstractionFlux
    changeStreamFlowCubicMetresPerSecond=self.perHourToPerSecond(changeStreamFlow*self.cellArea)
    return changeStreamFlowCubicMetresPerSecond

  def updateWaterHeightInChannel(self):
    self.waterHeightInChannel=(self.alpha*(self.Q**self.beta))/self.channelWidth;

  def updateWettedPerimeter(self):
    self.wettedPerimeter=wettedPerimeter(self.waterHeightInChannel,self.channelWidth)

  def updateAlpha(self):
    self.alpha=self.alpTerm*(self.wettedPerimeter**self.alpPow)
  
  def updateStore(self):
    self.store=(self.waterHeightInChannel*self.channelArea)/cellarea()

  def updateStates(self):
    self.updateWaterHeightInChannel()
    self.updateWettedPerimeter()
    self.updateAlpha()
    self.updateStore()

  def updateCumulativeDischarge(self):
    self.cumulativeDischargeCubicMetres=self.cumulativeDischargeCubicMetres + \
                                        self.perSecondToPerTimeStep(self.Q)

  def update(self,rainFlux,potentialAbstractionFlux):
    QIn=self.sideFlow(rainFlux,potentialAbstractionFlux)
    self.Q=kinematic(self.ldd,
                     self.Q,
                     QIn/self.distanceToDownstreamCell,
                     self.alpha,
                     self.beta,
                     ordinal(1),
                     self.timeStepDurationSeconds,
                     self.distanceToDownstreamCell)
    self.updateStates()
    self.updateCumulativeDischarge()
    return self.actualAbstractionFlux, self.Q

  def report(self, sample, timestep):
    report(self.Q,generateNameST('q', sample, timestep))
    report(self.waterHeightInChannel,generateNameST('h', sample, timestep))

  def reportInitial(self, sample):
    report(self.slopeToDownstreamNeighbour,generateNameS('std', sample))
    report(self.distanceToDownstreamCell,generateNameS('dtdc', sample))

