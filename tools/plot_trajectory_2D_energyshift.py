import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import glob
import sys

path = '/Users/tomkimpson/Data/ThesisData/RT/'
files = glob.glob(path + '*.txt')


plt.style.use('science')
fig, ax1,= plt.subplots(1, 1, figsize=(10,10))

def plot_it(f,id):

    data = np.loadtxt(f)
    x = data[:,0]
    y = data[:,1]
    nu = data[:,3]
    
    if id == 1:
        ax1.plot(x,y)

    else:
        t = np.linspace(0, 10, len(y))
        t = nu / nu[0]

        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = LineCollection(segments, cmap=plt.get_cmap('plasma'),
                            norm=plt.Normalize(0, 10))
        lc.set_array(t)
        plt.gca().add_collection(lc)
        ax1.scatter(x[0],y[0])


ID = 0
for f in files:
    plot_it(f,ID)

sq = 6
ax1.scatter(0,0,c='r')
plt.xlim(-sq, sq)
plt.ylim(-sq, sq)


fs = 20

#Label the axes
ax1.set_xlabel(r'$x [r_g]$',fontsize=fs)
ax1.set_ylabel(r'$y [r_g]$',fontsize=fs)

ax1.locator_params(axis='both', nbins=5) #set number of xticks                                                                                                                                                                           
ax1.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers






fname = 'RT_'+str(ID)+'.png'
savefile = '/Users/tomkimpson/Data/ThesisData/'+fname

dpi = 100
plt.savefig(savefile,dpi=dpi,bbox='tight')
plt.show()
