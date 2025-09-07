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

m=1.0
#N = 40
#vinit = 5.0
diameter = 1.0

def U(r,k):
    return 0.5*k * (diameter - r) **2

def Uavg(rhigh, rlow,k):
    #return k * (diameter - (rhigh - rlow)/2)
    return 1.0/6.0*k*((diameter-rhigh)**3.0-(diameter-rlow)**3.0)/(rlow-rhigh)

def rFromU(P,k):
    return diameter - math.sqrt(2.0*P/k) 

#def deltarstepping(deltar):
#    Nsteps = int(diameter / deltar)
#    return [diameter - i * deltar for i in range(Nsteps)]

def deltaUstepping(deltaU,k):
    rvals=[]
    for i in range(int(U(0,k) / deltaU)+1):
        r = rFromU(deltaU * i,k)
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

def average(rvals,k):
    Uvals=[0]
    for i in range(len(rvals)-1):
        rlow = rvals[i]
        rhigh = rvals[i+1]
        #print Uinte(rhigh, rlow)
    #################
        Uvals.append(Uavg(rhigh, rlow,k))
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
    evals=[]
    vel_values = []
    vel_times = []
    event_counter=0
    still_running = True
    while still_running:
        vijold = edmdsim.balls[0].v - edmdsim.balls[1].v
        vel_times.append(edmdsim.t)
        vel_values.append(numpy.linalg.norm(vijold))
        still_running = edmdsim.runEvent(vinit,re,deltaU,omega,beta,fraction)
        rij = edmdsim.balls[0].pos - edmdsim.balls[1].pos
        vij = edmdsim.balls[0].v - edmdsim.balls[1].v
        tcmax = math.pi / omega
        rcmax=vinit/omega*math.exp(-beta*tcmax/2.0)*math.sin(omega*tcmax/2.0)
        xvals.append(edmdsim.t/tcmax)            
        yvals.append((1.0-numpy.linalg.norm(rij))/rcmax)
        event_counter += 1

    return edmdsim.t,vij

print("Simulating tn")

def unconstantinftytc(deltaU,re,fraction,evals,k):
    global xvals, yvals
    rmax=1.0
#    n=int(0.5/deltaU)
    omega = math.sqrt(2 * k / m)
#fig = plt.figure()
#ax = fig.add_subplot(111)
#N = 1
#rmax=1.0
#deltar = rmax / N
#mu=1.0
#omega = math.sqrt(2 * k / m)
    gamma = -2 * math.log(re) * math.sqrt(k * 0.5 * m / (math.pi * math.pi + math.log(re) * math.log(re)))
#gamma=math.sqrt(2*k*m)
    beta = gamma / m
    edmdsim.lostenergy=0
    edmdsim.gamma = gamma
    omega02 = 2 * k / m
    if omega02 > beta * beta:
        omega = math.sqrt(omega02 - beta * beta)
    else:
        omega = math.sqrt(beta * beta - omega02)

    #tc = math.pi / omega
    #stepDiameter=deltarstepping(deltar=0.05)
    edmdsim.stepDiameter=deltaUstepping(deltaU,k)
    #print "step positions = ",edmdsim.stepDiameter
#stepEnergy=midSample(stepDiameter)
    edmdsim.stepEnergy=average(edmdsim.stepDiameter,k)
#    print "step energy = ",edmdsim.stepEnergy
    ShowPotential = False
    if ShowPotential:
        mp.step(stepDiameter, stepEnergy, where='pre')
        xvals = numpy.arange(0.0, 1.01, 0.01)
        mp.plot(xvals, U(xvals))
        mp.show()
    for e in evals:
        xvals=[]
        yvals=[]
        final=edmdsim.stepEnergy[len(edmdsim.stepEnergy)-1]
        unfinal=edmdsim.stepEnergy[len(edmdsim.stepEnergy)-4]
        KEinit=unfinal+e*(final-unfinal)
        vinit=2.0*math.sqrt(KEinit) 
        tmax,vij=collision_time(vinit,re,deltaU,omega,beta,fraction)

#        print "vinit",vinit
#        print "tmax",tmax
        #plt.plot(xvals,yvals, label=str(i))
        plt.plot(xvals,yvals, linestyle='-', marker='o', label=str(e), linewidth=2)
    #tcmax = 2.0*math.sqrt((stepEnergy[9]+stepEnergy[10])/2.0)*m/k
    #rcmax = -k/m*(0.5*tcmax)**2.0+2.0*math.sqrt((stepEnergy[9]+stepEnergy[10])/2.0)*0.5*tcmax
    #t = numpy.arange(0.0, 1.0, 0.01)
    #x=-4.0*t*t+4.0*t
    #plt.plot(t,x,'k--',linewidth=2)
    tcmax = math.pi / omega
    rcmax=vinit/omega*math.exp(-beta*tcmax/2.0)*math.sin(omega*tcmax/2.0)
    data = [(tv, vinit/omega*math.exp(-beta*tv)*math.sin(tv*omega)) for tv in numpy.arange(0.0, tcmax, 0.0001)]
    plt.plot([datum[0]/tcmax for datum in data], [datum[1]/rcmax for datum in data],'--',linewidth=2, color="black")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel('$t/t_c$',fontsize=18)
    plt.ylabel('$R/R_{max,c}$',fontsize=18)
    #plt.show()
#    print tmax
    return tmax,tcmax,vinit,KEinit

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

print("Fig 9")
unconstantinftytc(deltaU=0.01,re = 0.4,fraction=1.0/9.0,evals=[0.66, 0.665, 0.2, 0.9],k=1)
plt.annotate("$Steps\\approx22$",xy=(0.05,1.05), color="black")
plt.xlim(0,1.1)
plt.ylim(0,1.2)
plt.gcf().set_size_inches(5.0,3.6)
plt.tight_layout()
#plt.legend()
style_fig(plt.gcf())
plt.savefig("springdashpotSteps50e04traj.png", bbox_inches='tight',transparent=True, dpi=300)
plt.clf()

#print "Fig 10"
#unconstantinftytc(deltaU=0.01,re = 0.9,fraction=1.0/9.0,evals=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],k=1)
#plt.annotate("$Steps\\approx22$",xy=(0.025,1.0), color="black")
#plt.xlim(0,1.1)
#plt.ylim(0,1.1)
#plt.gcf().set_size_inches(5.0,3.6)
#plt.tight_layout()
##plt.legend()
#style_fig(plt.gcf())
#plt.savefig("springdashpotSteps50e09traj.png", bbox_inches='tight',transparent=True, dpi=300)
#plt.clf()

def KEindeltaU(deltaU,re,fraction,k,evals):
    global xvals, yvals
    rmax=1.0
#    n=int(0.5/deltaU)
    omega = math.sqrt(2 * k / m)
#fig = plt.figure()
#ax = fig.add_subplot(111)
#N = 1
#rmax=1.0
#deltar = rmax / N
#mu=1.0
#omega = math.sqrt(2 * k / m)
    gamma = -2 * math.log(re) * math.sqrt(k * 0.5 * m / (math.pi * math.pi + math.log(re) * math.log(re)))
#gamma=math.sqrt(2*k*m)
    beta = gamma / m
    edmdsim.lostenergy=0
    edmdsim.gamma = gamma
    omega02 = 2 * k / m
    if omega02 > beta * beta:
        omega = math.sqrt(omega02 - beta * beta)
    else:
        omega = math.sqrt(beta * beta - omega02)

    #tc = math.pi / omega
    #stepDiameter=deltarstepping(deltar=0.05)
    edmdsim.stepDiameter=deltaUstepping(deltaU,k)
    #print "step positions = ",edmdsim.stepDiameter
#stepEnergy=midSample(stepDiameter)
    edmdsim.stepEnergy=average(edmdsim.stepDiameter,k)
#    print "step energy = ",edmdsim.stepEnergy
    ShowPotential = False
    if ShowPotential:
        mp.step(stepDiameter, stepEnergy, where='pre')
        xvals = numpy.arange(0.0, 1.01, 0.01)
        mp.plot(xvals, U(xvals))
        mp.show()

    import numpy
    plot_e=[]
    plot_to=[]
    plot_ymax=[]
    plot_cor=[]
    for count, e in enumerate(evals):
        KEinit=e*deltaU
        vinit=2.0*math.sqrt(KEinit) 
        tcmax = math.pi / omega       
        tmax,vij=collision_time(vinit,re,deltaU,omega,beta,fraction)
        to=tmax/tcmax
        ro=max(yvals)
        vend=numpy.linalg.norm(vij)
        coeffre=vend/vinit
        plot_e.append(e)
        plot_to.append(to)
        plot_ymax.append(ro)
        plot_cor.append(coeffre)
        print('\r',e,count, int(count / float(len(evals)) * 10000) / 100.0,'%')

    plt.annotate("$e_{exact}="+str(re)+"$",xy=(35,1.2), color="black")
    plt.plot(plot_e,plot_to, linestyle='-', marker='', label=str(e), linewidth=1, color="black")
    
    #tcmax = 2.0*math.sqrt((stepEnergy[9]+stepEnergy[10])/2.0)*m/k
    #rcmax = -k/m*(0.5*tcmax)**2.0+2.0*math.sqrt((stepEnergy[9]+stepEnergy[10])/2.0)*0.5*tcmax
    #t = numpy.arange(0.0, 1.0, 0.01)
    #x=-4.0*t*t+4.0*t
    #plt.plot(t,x,'k--',linewidth=2)
#    plt.plot([datum[0]/tcmax for datum in data], [datum[1]/rcmax for datum in data],'--',color='black',linewidth=2)
    plt.plot([0,50],[1.0,1.0],'-', color='black')
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel('$KE_{init}/\\Delta{U}$ (Steps/Collision)',fontsize=18)
    plt.ylabel('$t_{max}/t_{max,c}$',fontsize=18)
    #plt.show()
#    print tmax
    plt.ylim(0.7,1.3)
    plt.gcf().set_size_inches(5.0,3.6)
    plt.tight_layout()
#plt.legend()
    style_fig(plt.gcf())
    plt.annotate("$e_{exact}="+str(re)+"$",xy=(35,1.2), color="black")
    plt.savefig("springdashpotSteps50e"+str(re)+"tc.png", bbox_inches='tight',transparent=True, dpi=300)
    plt.clf()
    
    plt.annotate("$e_{exact}="+str(re)+"$",xy=(35,1.2), color="black")
    plt.plot(plot_e, plot_ymax, linestyle='-', marker='', label=str(e), linewidth=1, color='black')
    tcmax = math.pi / omega
    rcmax=vinit/omega*math.exp(-beta*tcmax/2.0)*math.sin(omega*tcmax/2.0)
    plt.plot([0,50],[1.0,1.0],'-', color='black')
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel('$KE_{init}/\\Delta{U}$ (Steps/Collision)',fontsize=18)
    plt.ylabel('$R_{max}/R_{max,c}$',fontsize=18)
    #plt.show()
#    print tmax
    plt.ylim(0.7,1.3)
    plt.gcf().set_size_inches(5.0,3.6)
    plt.tight_layout()
#plt.legend()
    style_fig(plt.gcf())
    plt.savefig("springdashpotSteps50e"+str(re)+"rc.png", bbox_inches='tight',transparent=True, dpi=300)
    plt.clf()
            
    plt.plot(plot_e, plot_cor, linestyle='-', label=str(e), linewidth=1, color="black")
    plt.plot([0,50],[re,re],'-', color='black')
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel('$KE_{init}/\\Delta{U}$ (Steps/Collision)',fontsize=18)
    plt.ylabel('$e_{ED}$',fontsize=18)
    #plt.show()
#    print tmax
    plt.ylim(re-0.2,re+0.2)
    plt.gcf().set_size_inches(5.0,3.6)
    plt.tight_layout()
#plt.legend()
    style_fig(plt.gcf())
    plt.annotate("$e_{exact}="+str(re)+"$",xy=(35,re+0.1), color="black")
    plt.savefig("springdashpotSteps50e"+str(re)+"Cor.png", bbox_inches='tight',transparent=True, dpi=300)
    plt.clf()
    return tmax,tcmax,vinit,KEinit

print("Fig 11")
KEindeltaU(deltaU=0.01,re = 0.8,fraction=1.0/9.0,k=1,evals=[(i) for i in numpy.arange(1,50,0.005)])
plt.xlim(0,1.3)

