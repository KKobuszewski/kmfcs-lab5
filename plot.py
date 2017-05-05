from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mplcolors
import matplotlib.ticker



#Function to add ticks
def addticks(ax,newLocs,newLabels,pos='x',replace=True):
    # Draw to get ticks
    plt.draw()
    
    # Get existing ticks
    if pos=='x':
        locs   = ax.get_xticks().tolist()
        labels = [x.get_text() for x in ax.get_xticklabels()]
    elif pos =='y':
        locs = ax.get_yticks().tolist()
        labels=[x.get_text() for x in ax.get_yticklabels()]
    else:
        print("WRONG pos. Use 'x' or 'y'")
        return
    
    if replace is True:
        Dticks=dict(zip(newLocs,newLabels))
    else:
        # Build dictionary of ticks
        Dticks=dict(zip(locs,labels))
        
        # Add/Replace new ticks
        for Loc,Lab in zip(newLocs,newLabels):
            Dticks[Loc]=Lab

    # Get back tick lists
    locs=list(Dticks.keys())
    labels=list(Dticks.values())

    # Generate new ticks
    if pos=='x':
        ax.set_xticks(locs)
        ax.set_xticklabels(labels)
    elif pos =='y':
        ax.set_yticks(locs)
        ax.set_yticklabels(labels)



data = np.loadtxt('k-p.dat')

kx = data[:,0]
ky = data[:,1]
kz = data[:,2]
E1 = data[:,3]
E2 = data[:,4]
E3 = data[:,5]
E4 = data[:,6]

plt.title(r'$E\left( \mathbf{k} \right)$ for GaAs',fontsize=20)
plt.grid(True)
plt.xlabel(r'$k_z$',fontsize=15)
plt.ylabel('E',fontsize=15)

plt.plot(kz,E1,label='E1 (Ligth holes)',marker='o')
plt.plot(kz,E2,label='E2 (Ligth holes)',marker='.')
plt.plot(kz,E3,label='E3 (Heavy holes)',marker='x')
plt.plot(kz,E4,label='E4 (Heavy holes)',marker='*')

plt.legend(loc='lower left')
plt.savefig('kp.png')
plt.show()