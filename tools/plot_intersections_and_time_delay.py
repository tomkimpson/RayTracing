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


fig,(ax1,ax2)= plt.subplots(2, 1, figsize=(10,10),sharex=True)




#Load data
path = '/Users/tomkimpson/Data/ThesisData/RT/'


extra = 'grav_lensing_time_delay/'
primary_files = glob.glob(path +extra+ 'primary/*.txt')
secondary_files = glob.glob(path +extra+ 'secondary/*.txt')


#files = glob.glob(path + '*.txt')







def Format2D(ax):

    fs = 20
    #Label the axes

    ax.locator_params(axis='both', nbins=5)
    ax.tick_params(axis='both',which='major',labelsize=fs-4)









def plot(f,c):
    data = np.loadtxt(f)

    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    t = data[:,5]

    #Plot it
    #plot the rays
    ax1.scatter(x[-1],y[-1],c='0.1', marker = 'X')
    ax1.plot(x,y,c=c)

 #   and the time delays
    ax2.scatter(x[-1],t[-1],c=c)

    ax2.axvline(x[-1],linestyle='--',c='0.5')
    ax1.axvline(x[-1],linestyle='--',c='0.5')


for f in primary_files:
    plot(f,'C0')


for f in secondary_files:
    plot(f,'C2')

#The target point
#ax1.scatter(-150,60,c='C4', marker = 'X')


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20
Format2D(ax1)
Format2D(ax2)



ax2.set_xlabel(r'$x [r_g]$',fontsize=fs)
ax1.set_ylabel(r'$y [r_g]$',fontsize=fs)
ax2.set_ylabel(r'$t [r_g / c]$',fontsize=fs)
plt.setp(ax1.get_xticklabels(),visible=False)



#Draw the BH
ax1.scatter(0,0,c='r')


#axes limits
sq = 200
ax1.set_xlim(-180,10)
ax1.set_ylim(-10,180)

plt.subplots_adjust(hspace=0.01)

fname = 'grav_lensing_time_delay.png'
savefile = '/Users/tomkimpson/Data/ThesisData/'+fname
dpi = 100
plt.savefig(savefile,dpi=dpi,bbox='tight')


plt.show()



























