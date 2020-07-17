import csv, numpy

def frange(end,start=0,inc=0,precision=1):
    """A range function that accepts float increments."""
    import math

    if not start:
        start = end + 0.0
        end = 0.0
    else: end += 0.0

    if not inc:
        inc = 1.0
    count = int(math.ceil((start - end) / inc))

    L = [None] * count

    L[0] = end
    for i in (range(1,count)):
        L[i] = L[i-1] + inc
    return L

def getRecordFromCSVFile(csvFile,columnName,record):
  # reads a value from a table, columnNames are
  # given in the first row, first field in this row
  # has to be 'record'
  # record is the string in the first record column
  # returns a string
  fileReader=csv.DictReader(open(csvFile))
  records=[]
  rows=[]
  for row in fileReader:
    rows.append(row)
    records.append(row['record'])
  data=dict(zip(records,rows))
  return data[record][columnName]

def calculatePValueOfObservation(realizations,observation):
  # realizations: list with realizations of modelled distribution
  # observation: floating point value
  realizationsSorted=numpy.sort(realizations)
  #print 'sorted', realizationsSorted
  n=len(realizations)
  pValuesOfRealizations=(numpy.asarray(range(0,n,1))/float(n))+(0.5/n)
  #print 'pvalues', pValuesOfRealizations
  upperPosition=numpy.searchsorted(realizationsSorted,observation)
  lowerPosition=upperPosition-1
  if upperPosition == n:
    pValue=pValuesOfRealizations[n-1]
  elif upperPosition == 0:
    pValue=pValuesOfRealizations[0]
  else:
    pValue=pValuesOfRealizations[lowerPosition] + \
         (observation-realizationsSorted[lowerPosition])/(realizationsSorted[upperPosition]-realizationsSorted[lowerPosition]) *  \
         (pValuesOfRealizations[upperPosition]-pValuesOfRealizations[lowerPosition])
  return pValue


#realizations=[10,9,8,7,6]
##realizations=[0.2,0.9,1.0]
#observation=5.0
#p = calculatePValueOfObservation(realizations,observation)
#print p
