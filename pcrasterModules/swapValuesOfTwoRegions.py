import generalfunctions
from PCRaster import *

regions='mergeForestNoForest.map'
values='mergeVegAlbedoFS.map'

result=generalfunctions.swapValuesOfTwoRegions(regions,values,False)

report(result,'swapped.map')
