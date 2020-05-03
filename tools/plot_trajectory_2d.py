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




#Set up plotting environment
if (d == 3):
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1), (0,0),projection='3d')
elif  (d == 2):
    fig, ax1,= plt.subplots(1, 1, figsize=(10,10))




#Load data
path = '/Users/tomkimpson/Data/ThesisData/RT/'
files = glob.glob(path + '*.txt')







def Format2D(ax):

    fs = 20

    #Label the axes
    ax.set_xlabel(r'$x [r_g]$',fontsize=fs)
    ax.set_ylabel(r'$y [r_g]$',fontsize=fs)

    ax.locator_params(axis='both', nbins=5) #set number of xticks                                                                                                                                                                           
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers



    #Draw the BH
    ax.scatter(0,0,c='r')

    #axes limits
    sq = 300
    ax.set_xlim(-sq,sq)
    ax.set_ylim(-sq,sq)


def plot(f):
    data = np.loadtxt(f)

    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    #Plot it


    if (d == 3):
        ax1.plot(x,y,z)
        ax1.scatter(x[0],y[0],z[0], c='g')
        ax1.scatter(x[-1],y[-1],z[-1], c='r')
        ax1.scatter(0,0,0, c='k')


        limit = max(max(x),max(y),max(z))
        ax1.set_xlim(-limit,+limit)
        ax1.set_ylim(-limit,+limit)
        ax1.set_zlim(-limit,+limit)

    if (d == 2):
        ax1.plot(x,y)
        Format2D(ax1)


for f in files:
    plot(f)


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20



savefile = '/Users/tomkimpson/Data/ThesisData/raytracing_dispersion.png'
plt.savefig(savefile,dpi=300,bbox='tight')
plt.show()



























