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
fig, ax1,= plt.subplots(1, 1, figsize=(10,10))




#Load data
path = '/Users/tomkimpson/Data/ThesisData/RT/grav_lensing_time_delay/'

all_files = glob.glob(path + '../*.txt')
primary_files = glob.glob(path + 'primary/*.txt')
secondary_files = glob.glob(path + 'secondary/*.txt')



def Format2D(ax):

    fs = 20

    #Label the axes
    ax.set_xlabel(r'$\cos [r_g]$',fontsize=fs)
    ax.set_ylabel(r'$y [r_g]$',fontsize=fs)

    ax.locator_params(axis='both', nbins=5) #set number of xticks                                                                                                                                                                           
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers




    #axes limits
    sq = 200
    ax.set_xlim(-sq,sq)
 #   ax.set_ylim(-sq,sq)

def plot(f,c):
    data = np.loadtxt(f)

    phi = np.arccos(np.cos(data[-1,4]))
    
    x  = data[-1,0]
    t = data[-1,5]
        

    #Plot it
    ax1.scatter(x,t,c=c)


for f in all_files:
    plot(f,'C2')


#for f in primary_files:
 #   plot(f,'C0')

#for f in secondary_files:
 #   plot(f,'C2')

#The target point
#ax1.scatter(-150,60,c='C4', marker = 'X')


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20


Format2D(ax1)

fname = 'grav_lensing.png'
savefile = '/Users/tomkimpson/Data/ThesisData/'+fname
dpi = 100
#plt.savefig(savefile,dpi=dpi,bbox='tight')


plt.show()



























