#!/usr/bin/python

import edmdsim

import math 
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.pyplot as mp
import matplotlib.pyplot as plt
plt.rcParams['text.latex.preamble'] = "\\usepackage{amsmath}"
import numpy
k=1.0
m=1.0
diameter = 1.0

def U(r):
    return 0.5*k * (diameter - r) ** 2

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

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.clf()

def collision_time(vinit,re,deltaU,omega,beta,fraction):
    global capture_state,t, xvals, yvals, vel_values, vel_times    
    edmdsim.balls = [edmdsim.sphere(pos= (1.0, 0, 0), radius = 0.5 * diameter, m = m, v = edmdsim.vector(-0.5 * vinit, 0.0, 0.0)),
                     edmdsim.sphere(pos= (0.0, 0, 0), radius = 0.5 * diameter, m = m, v = edmdsim.vector(0.5 * vinit, 0.0, 0.0))]
    edmdsim.t=0.0
    capture_state = 0
    xvals =[]
    yvals =[]
    vel_values = []
    vel_times = []
    event_counter=0
    still_running = True
    while still_running:
        vijold = edmdsim.balls[0].v - edmdsim.balls[1].v
        vel_times.append(edmdsim.t)
        vel_values.append(vijold)
        still_running = edmdsim.runEvent(vinit,re,deltaU,omega,beta,fraction)
        rij = edmdsim.balls[0].pos - edmdsim.balls[1].pos
        vij = edmdsim.balls[0].v - edmdsim.balls[1].v
        tcmax = math.pi/math.sqrt(2*k/m)
        rcmax=vinit/math.sqrt(2*k/m)*math.sin(math.sqrt(2*k/m)*tcmax/2)
        xvals.append(edmdsim.t/tcmax)            
        yvals.append((1.0-numpy.linalg.norm(rij))/rcmax)
        event_counter += 1

    return edmdsim.t

def springonly(deltaU,re,fraction,evals):
    global xvals, yvals
    rmax=1.0
    n=int(0.5/deltaU)
    omega = math.sqrt(2 * k / m)
    gamma = -2 * math.log(re) * math.sqrt(k * 0.5 * m / (math.pi * math.pi + math.log(re) * math.log(re)))
    beta = gamma / m
    omega02 = 2 * k / m
    if omega02 > beta * beta:
        omega = math.sqrt(omega02 - beta * beta)
    else:
        omega = math.sqrt(beta * beta - omega02)
    edmdsim.stepDiameter=deltaUstepping(deltaU=deltaU)
    #print "step positions = ",edmdsim.stepDiameter
#stepEnergy=midSample(stepDiameter)
    edmdsim.stepEnergy=average(edmdsim.stepDiameter)
    #print "step energy = ",edmdsim.stepEnergy
    ShowPotential = False
    if ShowPotential:
        mp.step(stepDiameter, stepEnergy, where='pre')
        xvals = numpy.arange(0.0, 1.01, 0.01)
        mp.plot(xvals, U(xvals))
        mp.show()
    edmdsim.lostenergy=(1.0-re*re)*deltaU/(1.0+re*re)
    for e in evals:
        xvals=[]
        yvals=[]
        ten=edmdsim.stepEnergy[n]+edmdsim.lostenergy*(n)
        nine=edmdsim.stepEnergy[n-1]+edmdsim.lostenergy*(n-1)
        KEinit=nine+e*(ten-nine)
        vinit=2.0*math.sqrt(KEinit) 
        tmax=collision_time(vinit,re,deltaU,omega,beta,fraction)
        #plt.plot(xvals,yvals, label=str(i))
        plt.plot(xvals,yvals, linestyle='-', marker='o', label=str(e), linewidth=2)

    tcmax = math.pi/math.sqrt(2*k/m)
    rcmax=vinit/math.sqrt(2*k/m)*math.sin(math.sqrt(2*k/m)*tcmax/2)
    data=[(tv,vinit/math.sqrt(2*k/m)*math.sin(math.sqrt(2*k/m)*tv)) for tv in numpy.arange(0.0,tcmax,0.0001)]
    plt.plot([datum[0]/tcmax for datum in data], [datum[1]/rcmax for datum in data],'--',color='black',linewidth=2)
    plt.xlabel('$t/t_{c}$')
    plt.ylabel('$R/R_{max,c}$')
    return tmax,tcmax,vinit

def style_fig(figure):
    for axis in figure.axes:
        for i in ['bottom', 'top', 'left', 'right']:
            axis.spines[i].set_color('black')    
        axis.tick_params(axis='x', colors='black')
        axis.tick_params(axis='y', colors='black')
        axis.xaxis.label.set_color('black')
        axis.yaxis.label.set_color('black')
    factor=0.8
    figure.set_size_inches(5.0 * factor, 3.0 * factor)

springonly(deltaU=0.05,re=1.0,fraction=0.0,evals=[1-0.0001, 0+0.0001, 0.05,0.12,0.3])
plt.ylim(0,1.2)
plt.xlim(0,1.4)
plt.tight_layout()
plt.arrow(1.32, 1.01, 0.03, 0.0, fc='k', ec='k', head_length=0.03, facecolor=((0,0,0)), edgecolor=((0,0,0)))
plt.annotate("$\\infty$", xy=(1.32, 1.04), color="black")
plt.annotate("$\\tau=0^+$",xy=(1.0,1.01), color="black")
plt.annotate("$\\tau=0.05$",xy=(1.22,0.3),rotation=-65, color="black")
plt.annotate("$\\tau=0.12$",xy=(1.05,0.3),rotation=-65, color="black")  
plt.annotate("$\\tau=0.3$",xy=(0.9,0.3),rotation=-65, color="black")  
plt.annotate("$\\tau=1.0^-$",xy=(0.75,0.3),rotation=-65, color="black")
plt.annotate("$Steps\\approx8$",xy=(0.05,1.0), color="black")
#plt.legend()
style_fig(plt.gcf())
plt.savefig("spring_unfixed.png", bbox_inches='tight', dpi=300, transparent=True)
plt.clf()

springonly(deltaU=0.05,re=1.0,fraction=1.0/9.0,evals=[1.0/9-0.001,1.0/9+0.001, 0.3])
plt.ylim(0,1.2)
plt.xlim(0,1.2)
plt.tight_layout()
#    plt.arrow(1.32, 0.93, 0.03, 0.0, fc='k', ec='k', head_length=0.03)
plt.annotate("$\\tau=\\alpha^+$",xy=(1.05,0.3),rotation=-58, color='black')
plt.annotate("$\\tau=\\alpha^-$",xy=(0.72,0.3),rotation=-58, color='black')  
plt.annotate("$\\tau=0.3$",xy=(0.9,0.3),rotation=-58, color='black')
plt.annotate("$Steps\\approx8$",xy=(0.05,1.0), color="black")
    #plt.legend()
style_fig(plt.gcf())
plt.savefig("spring_fixed.png", bbox_inches='tight', dpi=300, transparent=True)
