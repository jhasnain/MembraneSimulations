import os 
import numpy as np
import matplotlib.pyplot as plt
import pymbar
import pickle
import scipy

bstr=4.0

data=open("tmp", "r").readlines()
dat=[]
for d in data:dat.append([float(x) for x in d.split()])
dat=np.array(dat)
bondnums = -dat[:, 8]/bstr

fit=[]
for i in range(5):fit.append(scipy.special.binom(4, i)*np.exp(bstr*i))

norm=(1+np.exp(bstr))**4
fit=np.array(fit)/norm

bins=np.array(range(-1, 5))+0.5
plt.clf()
plt.hist(bondnums, bins=bins, density=True)

nvals=np.array(range(0,5,1))
plt.plot(nvals, fit, marker='x')
plt.show()
