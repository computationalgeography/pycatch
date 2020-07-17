import numpy
import matplotlib.pyplot as plt

a=numpy.array([2,3,4])
print a
plt.plot(a,a+3,a,a+8,'g--')
plt.axis([2.5,3.5,2.2,18.3])
plt.show()
