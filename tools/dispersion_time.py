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


w = 20
h = 30



fig = plt.figure(figsize=(h,w))
ax1 = plt.subplot2grid((2,2), (0,0))
ax2 = plt.subplot2grid((2,2), (0,1),sharey=ax1)
ax3 = plt.subplot2grid((2,2), (1,0),colspan=2)


#Load data
path = '/Users/tomkimpson/Data/ThesisData/RT/dispersion/'
primary_rays = glob.glob(path + 'primary/*.txt')
secondary_rays = glob.glob(path + 'secondary/*.txt')






def plot(f,c):
    data = np.loadtxt(f)


    x = data[:,0]
    y = data[:,1]

    ax3.plot(x,y,c=c)

    
    t = data[:,5]
    nu = data[:,8]
    dt = t[-1] - t[0]

    code = data[0,7]


    return dt,nu[-1]/(2*np.pi*1e9),code









xx_primary=  []
yy_primary = []
for f in primary_rays:
    t,nu,code = plot(f,'C0')
    xx_primary.extend([t])
    yy_primary.extend([nu])


xx_secondary=  []
yy_secondary = []
for f in secondary_rays:
    t,nu,code = plot(f,'C2')
    xx_secondary.extend([t])
    yy_secondary.extend([nu])



#Normalise time
t0= min(xx_primary)
t1 = min(xx_secondary)


#Sort x,y arrays


Z = [x for _,x in sorted(zip(xx_primary,yy_primary))]
yy_primary = Z
xx_primary = sorted(xx_primary)


Z = [x for _,x in sorted(zip(xx_secondary,yy_secondary))]
yy_secondary = Z
xx_secondary = sorted(xx_secondary)











ax1.scatter(xx_primary-t0, yy_primary,c='C0')
ax1.scatter(xx_secondary-t1, yy_secondary,c='C2')

ax1.plot(xx_primary-t0, yy_primary,c='C0')
ax1.plot(xx_secondary-t1, yy_secondary,c='C2')



ax2.scatter(xx_primary-t0, yy_primary,c='C0')
ax2.scatter(xx_secondary-t0, yy_secondary,c='C2')

ax2.plot(xx_primary-t0, yy_primary,c='C0')
ax2.plot(xx_secondary-t0, yy_secondary,c='C2')


#X marks the spot
f = np.loadtxt(primary_rays[0])
xx = f[-1,0]
yy = f[-1,1]
ax3.scatter(xx,yy,c='0.1', marker = 'X')


#Prettyify
fs = 20

all_axes = plt.gcf().get_axes()
for ax in all_axes:
    ax.locator_params(axis='both', nbins=5) #set number of xticks
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers
    



#AX1 CONFIG
ax1.set_ylabel(r'$\nu$ [GHz]',fontsize=fs)
ax1.set_xlabel(r'$t [r_{\rm g}/c]$', fontsize=fs)
#AX2 CONFIG
plt.subplots_adjust(wspace=0.01)
plt.setp(ax2.get_yticklabels(),visible=False)
ax2.set_xlabel(r'$t [r_g/c]$',fontsize=fs)

#AX3 CONFIG

ax3.set_xlabel(r'$x [r_g]$', fontsize=fs)
ax3.set_ylabel(r'$y [r_g]$',fontsize=fs)
#axes limits
sq = 200
ax3.set_xlim(-200,400)
ax3.set_ylim(-50,150)
ax3.scatter(0,0,c='r')

#ax1.scatter(xx_primary - min(xx_primary), yy_primary)
#ax1.scatter(xx_secondary - t0_secondary, yy_secondary)


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20



fname = 'dispersion_time.png'
savefile = '/Users/tomkimpson/Data/ThesisData/'+fname
dpi = 300
plt.savefig(savefile,dpi=dpi,bbox='tight')


plt.show()



























