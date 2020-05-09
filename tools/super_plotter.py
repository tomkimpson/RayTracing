from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os
from matplotlib.collections import LineCollection



##Setup plotting environment
plt.style.use('science')

#Plot parameters
d = 2
plot_energy_shift = 'y' #'y', 'n' Only for 2D



#Set up plotting environment
if (d == 3):
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1), (0,0),projection='3d')
elif  (d == 2):
    fig, ax1,= plt.subplots(1, 1, figsize=(10,10))




#Load data
path = '/Users/tomkimpson/Data/ThesisData/RT/'
files = glob.glob(path + '*.txt')



def format_axes(ax):

    fs = 20

    #Label the axes
    ax.set_xlabel('x',fontsize=fs)
    ax.set_ylabel('y',fontsize=fs)


    #axes limits
    sq = 7
    

    ax.set_xlim(-sq,sq)
    ax.set_ylim(-sq,sq)

    #Just in 3D
    if d ==3:
        ax.set_zlabel('z',fontsize=fs)
        ax.set_zlim(-sq,sq)
        ax.view_init(elev=10., azim=148)
        ax.scatter(0,0,0,c='r')
    else:
        ax.scatter(0,0,c='r')

    ax.locator_params(axis='both', nbins=5)
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers


def plot(f):
    data = np.loadtxt(f)

    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    g_zamo = data[:,7] #g shift of zamo 
    Eobs = data[0,8]
    #Configure colors to convery gamma shift



    if (d == 3):
        ax1.plot(x,y,z)
        if (plot_energy_shift=='y'):
            print ('Note: you have selected to plot the energy shift in 3D mode which is currently not supported')


    if (d == 2):
        if plot_energy_shift=='n':

           ax1.plot(x,y)
        else:
            g_zamo = g_zamo/g_zamo[0]
            points = np.array([x, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, cmap=plt.get_cmap('plasma'),
                        norm=plt.Normalize(0, 10))
    
            lc.set_array(g_zamo)
            plt.gca().add_collection(lc)

#plot the rays
for f in files:
    plot(f)


format_axes(ax1) #pretty-ify
plt.show()



























