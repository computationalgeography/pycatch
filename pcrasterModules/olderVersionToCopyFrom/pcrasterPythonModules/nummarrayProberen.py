import numpy, math
#import scipy.linalg
#import scipy

#import pietje
#
##a = numpy.array([[4.0,0.1,0.5], \
##                 [0.1,2.0,0.7], \
##                 [0.5,0.7,2]])
##a = numpy.array([[4.0,1.1,1.5], \
##                 [1.1,2.0,1.7], \
##                 [1.5,1.7,2]])
#print 'covariantie matrix:'
#a = pietje.a
###########################################
# GEBRUIK numpy.loadtxt to read a array from disk.
# array with error of observations
# each row contains realizations of error of one observational variable
#data = numpy.array([[20, 10, 10.01, 20, 10],
#                    [10, 20, 10, 20.01,10],
#                    [10, 20, 10, 10, 20.01]
#                    ]
#                    )
#data = numpy.array([[10, 20, 10.01, 20, 10],
#                    [10, 20, 10, 20.01,10],
#                    [10, 20, 20, 10, 10.01]
#                    ]
#                    )
# covar werkt ook zo (1 variabele)
#data = numpy.array([[12,13,21,15,9]])

# conversions
# die covar heb ik al voor elke tijdstap
# functie nodig die die van disk leest
#test=numpy.corrcoef(data)
#print 'correlation matrix'
#print test
#print 'covar matrix'
#a = numpy.cov(data)
#print a
#print 'testje'
#a[0][0]=1
#a[1][1]=1
#a[2][2]=1

a=numpy.array(
   [[  1,   0.95,     0.95],
    [  0.95 ,  1 ,    0.95 ],
    [  0.95 ,  0.95 ,    1 ]]
      )

print 'covariace matrix'
print a
# inverse matrix
b = numpy.matrix(a).I
c = numpy.array(b)
print 'inverse matrix of covar matrix'
print c
# transpose matrix
#print a.transpose()
#print 'transpose matrix of covar matrix'
#print a.T

# observational data, length equals number of observational variables
# order of values should be equal to order of values in data (there: in rows)
# seems not to make a difference whether this is a column or row matrix
# --> functie nodig die uit set van kaarten (zelfde real/determ., t.s.) een array
# retourneert
#Hx = numpy.array([[2.0],[3.0],[5.0]])
#Hx = numpy.array([2.0,3.0,5.0])
#Hx = numpy.array([2.0])
# normale list mag ook
Hx = numpy.array([0.0,0.0,0.0])

# idem Hx, modelled data of realization
#y = numpy.array([[2.1],[3.5],[2.4]])
#y = numpy.array([2.1])
#y = numpy.array([4.0,1.0,2.0])

#################
# eerste particle
y = numpy.array([1 ,1 ,1])

# calculations
b=Hx-y
#print 'Hx'
#print Hx
#print 'y'
#print y
#print 'Hx-b'
#print b
#print 'transpose of (Hx-b):'
#print b.T
# firstTerm
firstTerm = numpy.dot(b.T,c)
#print 'firstTerm:'
#print firstTerm
#print 'wholeterm'
wholeTerm = numpy.dot(firstTerm,b)
#print wholeTerm
alphaFirstParticle= math.exp(-wholeTerm/2.0)

###############
# tweede particle
y = numpy.array([2,2,2])
# calculations
b=Hx-y
#print 'Hx'
#print Hx
#print 'y'
#print y
#print 'Hx-b'
#print b
#print 'transpose of (Hx-b):'
#print b.T
# firstTerm
firstTerm = numpy.dot(b.T,c)
print 'firstTerm:'
print firstTerm
print 'wholeterm'
wholeTerm = numpy.dot(firstTerm,b)
print wholeTerm
alphaSecondParticle=math.exp(-wholeTerm/2.0)

sumAlphas=alphaFirstParticle+alphaSecondParticle
weightFirstParticle=alphaFirstParticle/sumAlphas
weightSecondParticle=alphaSecondParticle/sumAlphas
print 'weight first particle is: ', weightFirstParticle
print 'weight second particle is: ', weightSecondParticle
print 'ratio of weights is: ', 1/(weightSecondParticle/weightFirstParticle)

