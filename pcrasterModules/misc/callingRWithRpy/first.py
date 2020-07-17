import rpy2.robjects as robjects
import numpy

#require=robjects.r['require']
#require("gstat")
##robjects.r('data(meuse)')
##a = robjects.r('''
##   variogram(zinc ~ 1, ~ x + y, meuse, boundaries=c(100,200,300,400))
##   ''')
##print a.names
##print a.r['dist']
##print type(a)
#
#z=range(1,10,1)
#zR=robjects.FloatVector(z)
#x=range(1,10,1)
#xR=robjects.FloatVector(x)
#y=range(1,10,1)
#yR=robjects.FloatVector(y)
#
#dataFrame=robjects.r['data.frame']
#framepje= dataFrame(x=xR,y=yR,z=zR)
#
#data=robjects.r['data']
#data((framepje))
#
#formulaOne=robjects.RFormula('z ~ 1')
#env=formulaOne.getenvironment()
#env['z']=zR
#formulaTwo=robjects.RFormula('~ x + y')
#env=formulaTwo.getenvironment()
#env['x']=xR
#env['y']=yR
#
##
##variogram=robjects.r['variogram']
##a=variogram(formulaOne, formulaTwo, framepje)
##
#
##print a.names

source=robjects.r['source']
source("test.r")
testFunction=robjects.r['testFunction']
a=testFunction(10)
print a
print type(a)

