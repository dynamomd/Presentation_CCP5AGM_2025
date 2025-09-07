#!/usr/bin/python

import math 
import numpy as np

def vector(a,b,c):
    return np.array([a,b,c])

balls = []
walls = []
stepDiameter=[]
stepEnergy=[]
NONE, INWARD, OUTWARD = range(3)
BallWallEVENT, BallBallEVENT, VirtualEVENT = range(3)
capture_state = 0
lostenergy=0
gamma=0
t = 0.0
infinity=float('inf')

def detectBallBallCollision(balli, ballj): 
    global capture_state
    rij = balli.pos - ballj.pos
    vij = balli.v - ballj.v
    rvdot = vij.dot(rij)
    diameter1 = stepDiameter[capture_state]
    c1 = rij.dot(rij) - diameter1 * diameter1 
    arg1 = rvdot * rvdot - vij.dot(vij) * c1
    if capture_state > 0:
        diameter2 = stepDiameter[capture_state-1]
        c2 = rij.dot(rij) - diameter2 * diameter2 
        arg2 = rvdot * rvdot - vij.dot(vij) * c2  
        if rvdot < 0 and arg1 > 0:
            return (-rvdot - math.sqrt(arg1))/(vij.dot(vij)), INWARD
        else:
            return (-rvdot + math.sqrt(arg2))/(vij.dot(vij)), OUTWARD
    elif capture_state == 0:
        if rvdot < 0 and arg1 > 0:            
            return (-rvdot - math.sqrt(arg1))/(vij.dot(vij)), INWARD
        else:
            return infinity, NONE
    else:
        raise Exception("Unknown capture_state")

def attractivecollision(balli, ballj, deltaKE, minEnergy):
    global capture_state
    rij = balli.pos - ballj.pos
    rij = rij / np.linalg.norm(rij)  
    vij = balli.v - ballj.v
    rvdot = vij.dot(rij)
    mu = balli.m * ballj.m / (balli.m + ballj.m)
    arg = rvdot * rvdot + 2 * deltaKE / mu
    
    if (arg - (deltaKE < 0) * 2 * minEnergy / mu) > 0:
        #capture in square well, association in square shoulder
        if rvdot < 0:            
            capture_state += 1
        elif rvdot > 0:
            capture_state -= 1
        else:
            raise Exception('Unexpected state')
        return 2 * rij * deltaKE / (rvdot+math.copysign(math.sqrt(arg), rvdot))
    else:
        ##bounce
        #tau = (0.5 * mu * rvdot * rvdot)/(stepEnergy[capture_state+1]-stepEnergy[capture_state])
        #laststep = capture_state
        #if tau > 1:
        #    tau = (0.5 * mu * rvdot * rvdot -(stepEnergy[capture_state+1]-stepEnergy[capture_state]))/(stepEnergy[capture_state+2]-stepEnergy[capture_state+1])
        #    laststep = capture_state+1
        #print "Bouncing with tau = ", tau, laststep
        return - rij * rvdot * 2 * mu

def runBallBallCollision(balli, ballj, eventtype, minenergy):
    global capture_state
    rij = balli.pos - ballj.pos
    rij = rij / np.linalg.norm(rij)  
    vij = balli.v - ballj.v
    rvdot = vij.dot(rij)

    if (capture_state == 0):
        deltar = 0
    else:
        deltar = stepDiameter[capture_state-1] - stepDiameter[capture_state]

    deltaKE = -lostenergy - gamma * abs(vij.dot(rij)) * deltar
    if eventtype == OUTWARD:
        if capture_state == 0:
            raise Exception('Went OUTWARD in step 0!')
        deltaKE += -(stepEnergy[capture_state - 1] - stepEnergy[capture_state])
        dP = attractivecollision(balli, ballj, deltaKE, minenergy)
        balli.v = balli.v + dP / balli.m
        ballj.v = ballj.v - dP / ballj.m
    else:
        if capture_state + 1 == len(stepEnergy):
            raise Exception("Heading INWARD on the innermost step!")
        deltaKE += -(stepEnergy[capture_state + 1] - stepEnergy[capture_state])
        dP = attractivecollision(balli, ballj, deltaKE, minenergy)
        balli.v = balli.v + dP / balli.m
        ballj.v = ballj.v - dP / ballj.m

def detectBallWallCollision(ball, wall):
    nw = wall.normal
    r = (ball.pos - wall.pos).dot(nw)
    if r < 0:
        nw = -nw
        r = -r
    v = (ball.v - wall.v).dot(nw)
    d = ball.radius + 0.5 * wall.size.dot(nw)
    if r < d and v < 0:
        return 0, 0
    if v >= 0:
        return infinity, 0
    return -(r - d) / v, 0

def runBallWallCollision(ball, wall, deltav):
    vij = ball.v - wall.v
    n = wall.normal
    ball.v = ball.v - 2 * vij.dot(n) * n

def stream(ball, deltat):
    ball.pos = ball.pos + ball.v * deltat

def streamSystem(deltat):
    global t
    t+= deltat
    for ball in balls:
        stream(ball, deltat)

def runEvent(vinit,re,deltaU,omega,beta,fraction):
    events=[]
    for balli in balls:
        for wall in walls:
            dt, deltav = detectBallWallCollision(balli, wall)
            events.append((dt, deltav, balli, wall, BallWallEVENT))
            
        for ballj in balls:
            if balli is not ballj:
                dt, deltav = detectBallBallCollision(balli, ballj)
                events.append(   (dt, deltav, balli, ballj, BallBallEVENT)    )
    
    nextevent = min(events)

    if nextevent[0] == infinity:
        return False #This means we're out of events!

    streamSystem(nextevent[0])

    if nextevent[4] == BallWallEVENT:
        runBallWallCollision(nextevent[2], nextevent[3], nextevent[1])
    elif nextevent[4] == BallBallEVENT:
        runBallBallCollision(nextevent[2], nextevent[3], nextevent[1], fraction * deltaU)
    else:
        raise Exception("Unknown event type")
    return True

class sphere:
    def __init__(self, pos, radius, m, v):
        self.pos = np.array(pos)
        self.radius = radius
        self.m = m
        self.v = np.array(v)  
    
    def __lt__(self, other):
        return self.pos[0] < other.pos[0]

def calc_rel_ke_energy(balli,ballj):
    sumKE = 0.0
    vij=balli.v-ballj.v
    for ball in balls:
        sumKE += vij.dot(vij) * 0.5 * ball.m    
    return sumKE, stepEnergy[capture_state]

def calc_energy():
    sumKE = 0.0
    for ball in balls:
        sumKE += ball.v.dot(ball.v) * 0.5 * ball.m
    return sumKE, stepEnergy[capture_state]