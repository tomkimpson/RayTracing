from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os



##Setup plotting environment
plt.style.use('science')
d = 2

fig, ax1,= plt.subplots(1, 1, figsize=(10,10))




#Load data
path = '/Users/tomkimpson/Data/ThesisData/RT/'
f = path + 'optimisation.txt'



def Format2D(ax):

    fs = 20

    #Label the axes
    ax.set_xlabel(r'$x [r_g]$',fontsize=fs)
    ax.set_ylabel(r'$y [r_g]$',fontsize=fs)

 #   ax.locator_params(axis='both', nbins=5) #set number of xticks                                                                                                                                                                           
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers



    #Draw the BH

def plot(f):
    data = np.loadtxt(f)

    x = data[:,0]
    y = data[:,1]

    ax1.plot(x,y)
    Format2D(ax1)


plot(f)


#The target point


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20


plt.show()



























