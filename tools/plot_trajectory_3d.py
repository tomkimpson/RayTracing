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
d = 3




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
    ax.set_xlabel('x',fontsize=fs)
    ax.set_ylabel('y',fontsize=fs)

    ax.locator_params(axis='both', nbins=5) #set number of xticks                                                                                                                                                                           
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers



    #Draw the BH
    ax.scatter(0,0,c='r')

    #axes limits
    sq = 7
    ax.set_xlim(-sq,sq)
    ax.set_ylim(-sq,sq)



def Format3D(ax):

    fs = 20

    #Label the axes
    ax.set_xlabel('x',fontsize=fs)
    ax.set_ylabel('y',fontsize=fs)
    ax.set_zlabel('z',fontsize=fs)

    ax.locator_params(axis='both', nbins=5) #set number of xticks                                                                                                                                                                           
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers



    #Draw the BH
    ax.scatter(0,0,0,c='r')

    #axes limits
    sq = 4
    ax.set_xlim(-sq,sq)
    ax.set_ylim(-sq,sq)
    ax.set_zlim(-sq,sq)


    ax.view_init(elev=10., azim=148)

#    ax.grid(False)
#    ax.xaxis.pane.set_edgecolor('black')
#    ax.yaxis.pane.set_edgecolor('black')
#    ax.xaxis.pane.fill = False
#    ax.yaxis.pane.fill = False
#    ax.zaxis.pane.fill = False



def plot(f):
    data = np.loadtxt(f)

    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    #Plot it


    if (d == 3):
        ax1.plot(x,y,z)

    if (d == 2):
        ax1.plot(y,z)


for f in files:
    plot(f)






if (d==3):
    Format3D(ax1)
if (d==2):
    Format2D(ax1)







plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20



savefile = '/Users/tomkimpson/Data/ThesisData/raytracing3d.png'
plt.savefig(savefile,dpi=1200,bbox='tight')

plt.show()



























