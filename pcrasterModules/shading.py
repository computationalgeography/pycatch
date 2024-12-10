import sys
import pysolar.solar
import pytz
import math
import supportingfunctions
from pcraster.framework import *
from pcraster import *
from datetime import *
import datetimePCRasterPython
import component

# notes
# time step duration in h
# vertical fluxes in m/h, variable name 'flux'
# vertical fluxes over a time step in m, variable name 'fluxAmount'
# amounts in storages in m (stores), variable name 'store'
# (everything as a waterslice over a whole cell)
# except where indicated
# inputs of functions may be python types, return values of
# functions are always PCRaster types

# note that south should be bottom of map!!
# note that units of x coor should be equal to units of elevation on digital elevation map

def createListOfSolarCritAngles(step, extendedDem):
    # creates a list where each element is a list with two items: first
    # item is the azimuth, second item is the solar crit angle map for
    # that azimuth
    # step is a floating point, recommended e.g. 0.1, ie about 6 degrees
    azimuthRadians=supportingfunctions.frange(2*math.pi, 0.0, step, 2) + [2*math.pi]
    set=[]
    for azimuth in azimuthRadians:
        # horizontan, assumed here it returns angle in radians, 1.57 is 90 degrees, position of sun
        # 90 degrees; the cover operation is needed as horizontan does not calculate angles on the edge
        # of the map, it is assumed these can always be in the sun (during the day)
        solarCritAngle=cover(horizontan(extendedDem, azimuth), 1.57)
        set.append([azimuth, solarCritAngle])
    return set

class Shading(component.Component):
    def __init__(self, digitalElevationModel, latitudeFloatingPoint, longitudeFloatingPoint, timeZone, reduceRunTime, \
                timeStepsToReport, setOfVariablesToReport):

        # init for supsend and resume in filtering only
        self.fractionSolarBeamReceived=scalar(0)
        self.shaded=scalar(0)
        self.shadedNoHorizonEffect=scalar(0)
        self.solarAltitudeDegrees=scalar(0)
        self.solarAltitudeRadians=scalar(0)
        self.solarAzimuthDegrees=scalar(0)
        self.solarAzimuthRadians=scalar(0)
        self.solarAzimuthRadiansConverted=scalar(0)
        self.solarCritAngle=scalar(0)
        self.variablesToReport={}

        # real inits

        # dem: bottom should be south
        # reduceRunTime results in precalculation of horizontan output, with a step of 0.1 (see
        # createListOfSolarCritAngles), the difference is really small because it only has effect on the
        # shading and this is a very small number of pixels and mostly with a low solar angle already anyway
        self.digitalElevationModel=digitalElevationModel
        extendedDem=windowaverage(self.digitalElevationModel, celllength()*3)
        self.extendedDem=cover(self.digitalElevationModel, extendedDem)
        self.latitudeFloatingPoint=latitudeFloatingPoint
        self.longitudeFloatingPoint=longitudeFloatingPoint
        self.timeZone=timeZone
        self.demSlopeAngle=atan(slope(self.digitalElevationModel))
        self.aspect=aspect(digitalElevationModel)
        self.aspectWithFlat=ifthenelse(nodirection(self.aspect), 1.0, self.aspect)
        self.reduceRunTime=reduceRunTime
        if self.reduceRunTime:
            self.solarCritAngleList=createListOfSolarCritAngles(0.1, self.extendedDem)

        self.timeStepsToReport=timeStepsToReport
        self.setOfVariablesToReport=setOfVariablesToReport

    def reportAsMaps(self, sample, timestep):
        self.output_mapping = {
                               'Mfs': self.fractionSolarBeamReceived,
                               'Msc': self.solarCritAngle,
                               'Msh': self.shaded
                              }
        self.variablesToReport = self.rasters_to_report(self.setOfVariablesToReport)
        self.reportMaps(sample, timestep)

    def updateVariablesAsNumpyToReport(self):
        self.variablesAsNumpyToReport = {
                                        }

    def reportAsNumpyOneFile(self, locations, sample, timestep, endTimeStep):
        self.updateVariablesAsNumpyToReport()
        self.reportAsNumpyOneFilePerRealization(locations, sample, timestep, endTimeStep)

    def reportAsNumpyMultipleFiles(self, locations, sample, timestep):
        self.updateVariablesAsNumpyToReport()
        self.reportAsNumpy(locations, sample, timestep)

    def horizontanFromList(self, azimuth):
        # retrieve critical solar angle from list
        for item in self.solarCritAngleList:
            if item[0] > azimuth:
                result = item
                break
        solarCritAngle=result[1]
        return solarCritAngle

    def calculateSolarAngles(self, timeDatetimeFormat):
        # convert naive time to aware time
        timezone = pytz.timezone(self.timeZone)
        d_aware = timezone.localize(timeDatetimeFormat)
        print(timeDatetimeFormat, d_aware)
        # solar altitude, sometimes called solar angle
        self.solarAltitudeDegrees=pysolar.solar.get_altitude(self.latitudeFloatingPoint, self.longitudeFloatingPoint, d_aware)
        self.solarAltitudeRadians=math.radians(self.solarAltitudeDegrees)
        # azimuth: west is -90 degrees, north is -180 degrees, east is -270 degrees
        self.solarAzimuthDegrees=pysolar.solar.get_azimuth(self.latitudeFloatingPoint, self.longitudeFloatingPoint, d_aware)
        self.solarAzimuthRadians=0.0-math.radians(self.solarAzimuthDegrees)
        # solarAzimuthRadiansConverted: north is 0, east is 0.5 pi, south is pi, west is 1.5 pi
        if self.solarAzimuthRadians < math.pi:
            self.solarAzimuthRadiansConverted=self.solarAzimuthRadians+math.pi
        else:
            self.solarAzimuthRadiansConverted=self.solarAzimuthRadians-math.pi
        if self.reduceRunTime:
            self.solarCritAngle=self.horizontanFromList(self.solarAzimuthRadiansConverted)
        else:
            self.solarCritAngle=horizontan(self.extendedDem, self.solarAzimuthRadiansConverted)

    def updateShading(self, timeDatetimeFormat):
        self.shadedNoHorizonEffect=pcrgt(scalar(self.solarCritAngle), math.radians(self.solarAltitudeDegrees))
        self.shaded=ifthenelse(pcrlt(scalar(self.solarAltitudeDegrees), scalar(0.0)), boolean(1), self.shadedNoHorizonEffect)
        return self.solarCritAngle, self.shaded

    def updateIncidence(self, timeDatetimeFormat):
        # incidence is the angle between the perpendicular plane of the incoming solar rays and the
        # surface on which they are projected (i.e. the digital elevation model)
        # equation 5 in O. van Dam, thesis, p. 71, note that last term should read cos(X) instead of cos(B)
        termOne=math.cos(self.solarAltitudeRadians)*sin(self.demSlopeAngle)
        a=self.solarAzimuthRadiansConverted
        b=self.aspectWithFlat
        termTwo=scalar(math.sin(a))*sin(b)+scalar(math.cos(a))*cos(b)
        termThree=math.sin(self.solarAltitudeRadians)*cos(self.demSlopeAngle)
        self.fractionSolarBeamReceived=max(termOne * termTwo + termThree , 0.0)
        fractionSolarBeamReceivedFlatSurfaceAlsoNegative=math.sin(self.solarAltitudeRadians)
        # if statement required because max([a,b]) werkt niet na importeren PCRaster
        if fractionSolarBeamReceivedFlatSurfaceAlsoNegative < 0.0:
            fractionSolarBeamReceivedFlatSurface=0.0
        else:
            fractionSolarBeamReceivedFlatSurface=fractionSolarBeamReceivedFlatSurfaceAlsoNegative
        return self.fractionSolarBeamReceived, fractionSolarBeamReceivedFlatSurface

    def update(self, timeDatetimeFormat):
        self.calculateSolarAngles(timeDatetimeFormat)
        self.solarCritAngle, self.shaded=self.updateShading(timeDatetimeFormat)
        fractionSolarBeamReceivedNoShading, fractionSolarBeamReceivedFlatSurface =self.updateIncidence(timeDatetimeFormat)
        fractionSolarBeamReceivedWithShading=ifthenelse(self.shaded, scalar(0), fractionSolarBeamReceivedNoShading)
        return fractionSolarBeamReceivedWithShading, fractionSolarBeamReceivedFlatSurface, self.shaded

    def printit(self):
        print('solar altitude: ', self.solarAltitudeDegrees, 'solar azimuth: ', self.solarAzimuthDegrees)
        print('solar azimuth radians, South is pi: ', self.solarAzimuthRadiansConverted)

## test
#
## all time UTC
#
#class ShadingModel(object):
#  def __init__(self):
#    pass
#
#  def premcloop(self):
#    self.d_dateTimePCRasterPython=datetimePCRasterPython.DatetimePCRasterPython(datetime(2009,7,27),timedelta(0,0,0,0,0,0.1,0))
#    self.d_shading=Shading(scalar('mdtpaz4.map'),52.1283333333333 , 5.19861111111111)
#
#  def initial(self):
#    self.totalFraction=scalar(0)
#
#  def dynamic(self):
#    print self.currentTimeStep()
#    self.d_dateTimePCRasterPython.update()
#    self.d_dateTimePCRasterPython.printit()
#    timeDatetimeFormat=self.d_dateTimePCRasterPython.getTimeDatetimeFormat()
#    fractionReceived=self.d_shading.update(timeDatetimeFormat)
#    self.totalFraction=self.totalFraction+fractionReceived
#    self.d_shading.printit()
#    self.d_shading.report(self.currentSampleNumber(), self.currentTimeStep())
#
#  def postmcloop(self):
#    pass
#    self.report(self.totalFraction,'total')
#
#myModel = ShadingModel()
#dynamicModel = DynamicFramework(myModel, 30*24)
#mcModel = MonteCarloFramework(dynamicModel, 1)
#mcModel.run()
