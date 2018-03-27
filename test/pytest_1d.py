#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 11:50:25 2018

@author: votec
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn import gaussian_process
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
import time
import sklearn
import gp

def obj_fun(xs):
#    return np.sin(xs) + np.random.randn(1)*0.05
    return np.exp(-(xs - 2)**2) + np.exp(-(xs - 6)**2/10) + 1/ (xs**2 + 1)

def rosenbrock(x, dim):
    y = 0
    for i in range(dim-1):
        y += 100 * (x[i+1] - x[i] * x[i])**2 + (1 - x[i])**2
    return y

x = np.linspace(-2, 10, 100)
y = obj_fun(x)
plt.plot(x, y)

x0 = np.linspace(-2, 10, 130).reshape(-1,1)
y0 = obj_fun(x0)
#plt.plot(x0, y0, 'x')
#kernel = GPy.kern.RBF(input_dim=1, variance=1., lengthscale=1.)
s = time.time()
gpr = GaussianProcessRegressor(
            kernel=Matern(),
#            kernel=sklearn.gaussian_process.kernels.RBF(),
#            optimizer=None,
            normalize_y = False,
            n_restarts_optimizer=0
        )
#gp = gaussian_process.GaussianProcess(theta0=1e-2, thetaL=1e-4, thetaU=1e-1)
gpr.fit(x0, y0)

x1 = np.linspace(-2, 10, 100).reshape(-1,1)
y1, std = gpr.predict(x.reshape(-1,1),return_std=True)
y1 = y1.reshape(-1)
std= std.reshape(-1)
plt.plot(x, y1, 'r-')
plt.fill_between(x, y1-std, y1+std, alpha=0.2)
print(time.time() - s)
print('mean', np.mean(np.abs(y1-y)))
print('max',  np.max(np.abs(y1-y)))
print('argmax',  x[np.argmax(np.abs(y1-y))],y[np.argmax(np.abs(y1-y))],y1[np.argmax(np.abs(y1-y))])



s = time.time()
a=gp.gp(order=1, norm=True)
a.fit(x0, y0)
yc, std=a.predict(x.reshape(-1,1), return_std=True)
yc, std = np.array(yc).reshape(-1),np.array(std).reshape(-1)
print(time.time() - s)
plt.plot(x, yc, 'g-')
plt.fill_between(x, yc-std, yc+std, alpha=0.2)

#plt.plot(x, yc[1], '-')
plt.show()

