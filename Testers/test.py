import numpy as np


x = np.array([0, 1, 2, 3])
y = 2*x + 1

print y
sum = 0

for i in xrange(x.size-1):
    sum += (x[i+1]-x[i])*(y[i]+y[i+1])/2

print sum
