#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 11:50:25 2018

@author: votec
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import gaussian_process
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
import time
import sklearn
import gp

DIM = 2
N_SAMPLE = 100
TEST_SAMPLE = 50 # ^2

def rosenbrock(x, dim):
    y = 0
    for i in range(dim-1):
        y += 100 * (x[i+1] - x[i] * x[i])**2 + (1 - x[i])**2
    return y

def show(y):
    plt.imshow(np.log(y).reshape(TEST_SAMPLE,TEST_SAMPLE),extent=(-2,2,-2,2))
    plt.show()
def show3d(y, x):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(x[:,0].reshape(TEST_SAMPLE,TEST_SAMPLE), x[:,1].reshape(TEST_SAMPLE,TEST_SAMPLE), y.reshape(TEST_SAMPLE,TEST_SAMPLE), rstride=1, cstride=1, cmap='rainbow')
    ax.scatter(x0[:,0].reshape(-1), 
               x0[:,1].reshape(-1), 
               y0.reshape(-1), c='r')
    ax.view_init(azim=60)
    plt.show()
x0 = np.random.random((N_SAMPLE, DIM)) * 4 - 2
y0 = rosenbrock(x0.T, DIM).reshape(-1,1)
x = np.linspace(-2, 2, TEST_SAMPLE)
x = np.array(np.meshgrid(x,x)).T.reshape(TEST_SAMPLE**2, 2)
y = rosenbrock(x.T, DIM)
#show3d(y,x)
gpr = GaussianProcessRegressor(
            kernel=Matern(),
#            kernel=sklearn.gaussian_process.kernels.RBF(),
            optimizer=None,
            normalize_y = False,
            n_restarts_optimizer=0
        )
gpr.fit(x0, y0)
ypy = gpr.predict(x)
#show(y_t)
show3d(ypy,x)
gpr = GaussianProcessRegressor(
            kernel=Matern(),
#            kernel=sklearn.gaussian_process.kernels.RBF(),
#            optimizer=None,
            normalize_y = False,
            n_restarts_optimizer=0
        )
gpr.fit(x0, y0)
ypy = gpr.predict(x)
#show(y_t)
show3d(ypy,x)
a=gp.gp(order=1, norm=True)
a.fit(x0, y0)
yc = a.predict(x)
yc = np.array(yc[0])
show3d(yc,x)
a=gp.gp(order=0, norm=True)
a.fit(x0, y0)
yc = a.predict(x)
yc = np.array(yc[0])
show3d(yc,x)
#yc, std=a.predict(x, return_std=True)