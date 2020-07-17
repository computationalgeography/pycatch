import math, supportingfunctions

azimuthRadians=supportingfunctions.frange(2*math.pi,0.0,0.1,2)

set=[]
i=1
for azimuth in azimuthRadians:
  set.append([azimuth,i])
  i = i+1

print set

value=0.0
for item in set:
  if item[0] > value:
    result = item
    break

secondValue=result[1]
print secondValue
    

  
