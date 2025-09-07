#!/usr/bin/python
import random
import pylab
from matplotlib import rc
import numpy as np
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

import matplotlib.pyplot as mp
import matplotlib.pyplot as plt
import math
import numpy
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['text.latex.preamble'] = "\\usepackage{amsmath}"

#########MAIN PLOT

k=1.0
m=1.0
N = 10
rmax=1.0
omega = math.sqrt(2 * k / m)
#N = 1
#rmax=1.0
#deltar = rmax / N
k=1.0
#vinit = 5.0
diameter = 1.0
m = 1.0
deltaU=0.05
#mu=1.0
#omega = math.sqrt(2 * k / m)

def U(r):
    return 0.5*k * (diameter - r)**2.0 

def Uavg(rhigh, rlow):
    #return k * (diameter - (rhigh - rlow)/2)
    return 1.0/6.0*k*((diameter-rhigh)**3.0-(diameter-rlow)**3.0)/(rlow-rhigh)

def rFromU(P):
    return diameter - math.sqrt(2.0*P/k)

#def deltarstepping(deltar):
#    Nsteps = int(diameter / deltar)
#    return [diameter - i * deltar for i in range(Nsteps)]

def deltaUstepping(deltaU):
    rvals=[]
    for i in range(int(U(0) / deltaU)+1):
        r = rFromU(deltaU * i)
        if r != 0:
            rvals.append(r)
  
    return rvals+[0]
    
#def midSample(rvals):
#    Uvals=[0]
#    for i in range(len(rvals)-1):
#        rlow = rvals[i]
#        rhigh = rvals[i+1]
#        Uvals.append(U((rlow+rhigh)/2.0))
#    return Uvals

def average(rvals):
    Uvals=[0]
    for i in range(len(rvals)-1):
        rlow = rvals[i]
        rhigh = rvals[i+1]
        #print Uinte(rhigh, rlow)
    #################
        Uvals.append(Uavg(rhigh, rlow))
    #################
    return Uvals

#stepDiameter=deltarstepping(deltar=0.05)
stepDiameter=deltaUstepping(deltaU=deltaU)
print("step positions = ",stepDiameter)
#stepEnergy=midSample(stepDiameter)
stepEnergy=average(stepDiameter)
print("step energy = ",stepEnergy)
#stepDiameter=[1,0.5,0.5]
#stepEnergy=[0,0,1]

plt.step(stepDiameter, stepEnergy, where='pre',c='green',linewidth=2)
xvals = numpy.arange(0.0, 1.05, 0.01)
plt.plot(xvals, U(xvals),c='blue',linewidth=2)
plt.plot([0.0,1.1],[0.35,0.35],'--',c='black',linewidth=2)
#plt.plot([1.0,1.0,1.0],[0.0,0.05,0.1],'--',c='c',linewidth=2)
plt.xticks([0.0,0.2,0.4,0.6,0.8,1.0],['','','','','','$R_{cutoff}$'])
plt.yticks([0.0,0.1,0.2,0.3,0.35,0.4,0.5],['','','','','$KE_{in}$','',''])
plt.ylabel('$U_{i}$')
plt.xlabel('$R$')
plt.xlim(0,1.05)

def style_fig(figure):
    for axis in figure.axes:
        for i in ['bottom', 'top', 'left', 'right']:
            axis.spines[i].set_color('black')    
        axis.tick_params(axis='x', colors='black')
        axis.tick_params(axis='y', colors='black')
        axis.xaxis.label.set_color('black')
        axis.yaxis.label.set_color('black')

style_fig(plt.gcf())
a = plt.axes([0.45, 0.45, 0.4, 0.4])
style_fig(plt.gcf())


################SUBPLOT
def yc(x):
    return 1.2-x
xx=[1.2, 0.4, 0.4, 0.0]
yy=[0.0, 0.4, 1.2, 1.2]
datac=[(x,yc(x)) for x in numpy.arange(0.0,1.2,0.001)]

#plt.clf()
plt.plot([data[0] for data in datac], [data[1] for data in datac],linewidth=2, c="black")
plt.step(xx,yy,linewidth=2)
#plt.plot([0.4,0.4,0.4],[0.0,0.2,0.4],'--')
plt.plot([0.0, 0.2, 0.4],[0.4, 0.4, 0.4], '--', c='black', linewidth=3)
#plt.plot([0.0, 0.6, 1.2],[0.6, 0.6, 0.6], '--',linewidth=4)
#plt.plot([0.9, 0.9, 0.9],[0.4, 0.5, 0.6], '--')
plt.plot([0.4,0.4,0.4],[0.0,0.2,0.4],'--', c='black', linewidth=2)
plt.xlim(0.0, 1.4)
plt.ylim(0.0, 1.4)
labely = [item.get_text() for item in a.get_yticklabels()]
labelx=[item.get_text() for item in a.get_xticklabels()]
#labelx[5]='$\lambda\sigma$'
#labelx[7]='$R$'
#labelx[2]='$r_{final}$'
labelx[2]='$R_{final+1}$'
#labelx[6]='$r_{final+1}$'
labelx[6]='$R_{final}$'
labely[2]='$U_{final}$'

labely[6]='$U_{final+1}$'
factor=0.8
plt.gcf().set_size_inches(5.0 * factor, 3.0 * factor)
#plt.annotate('$KE_{final}$', (0.9,0.6), xytext=(0.9,0.4), arrowprops=dict(facecolor='red'
#	       ),fontsize=17, horizontalalignment='left', verticalalignment='top')
plt.annotate("", xy=(0.85, 0.38), xytext=(0.85, 0.7), xycoords="data", textcoords="data", arrowprops=dict(arrowstyle='<->', facecolor=((0,0,0)), edgecolor=((0,0,0))), color="black")
plt.annotate("$KE_{final}$", xy=(0.9, 0.45), color="black")
#ax.set_xlabel('$R$', fontsize = 13)
#ax.set_ylabel('$U_{i}$', fontsize = 13)
#plt.tight_layout()
a.set_yticklabels(labely)
a.set_xticklabels(labelx)
plt.savefig("single_sim.png", pad_inches=0, transparent=True, bbox_inches='tight', dpi=300)

exit()
