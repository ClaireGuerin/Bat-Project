# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 14:31:31 2016

@author: Claire
"""

import random as rd
import math  as m
import numpy as np
import matplotlib.pyplot as plt

###----------LAUNCHER----------###
### Class implementation to launch the simulation environment.
### By convention, an object of class Launcher() is stored under the name env (stands for environment)
### Initiating inputs: popsize, boxsize, simduration, realduration
### popsize: positive integer. number of agents to simulate i.e. size of the bat group.
### boxsize: list of 2 elements. Cartesian space within which the agents can move
### simduration: positive integer. Duration of simulation 
### ### sim stands for simulation throughout code
### realduration: positive integer. "Real" duration the simulation is supposed to reflect. 

###----------INITIATING OUTPUTS----------###
### all_ID: empty array of dimensions 1x(n=popsize elements). 
### all_initpos: empty list of n=popsize lists of 2 elements each. 
### ### init & pos respectively stand for initial and position throughout the code
### timecount: positive integer. Timing of simulation iteration.
### timetracking: list of 1 element. Initial timing of the simulation
### timeresolution: positive integer. Resolution of the time simulated.

#----------IDENTIFICATION----------#
# Function listing the agents by ID, from 0 to popsize-1.
# Inputs: popsize taken from __init__ method.

#----------OUTPUTS----------#
# all_ID: array of dimensions 1 x (n=popsize). ID of each agent from 0 to n-1.
#----------end of documentation for IDENTIFICATION----------#

#----------POSITIONS--------#
# Function randomly attributing initial positions to each agent
# Inputs: all_ID taken from Identification method.

#----------OUPTPUTS----------#
# all_initpos: updated all_initpos. Contains initial position(s) of agent(s)
#----------end of documentation for POSITIONS----------#

#----------TIMELINE----------#
# Function drawing the timeline of the simulation
# Inputs: taken from __init__ method

#----------OUTPUTS----------#
# timecount: updated timecount. At the end of the simulation, timecount = simduration / timeresolution
# timetracking: array of dimensions 1xsimduration. Updated timetracking with "Real times" corresponding to each simulation iteration
#----------end of documentation for TIMELINE----------#
###----------end of documentation for LAUNCHER----------#  

class Launcher:
    
    def __init__(self, popsize, boxsize, simduration, realduration):
        self.popsize = popsize
        self.simduration = simduration
        self.timecount = 0
        self.timetracking = [self.timecount]
        self.realduration = realduration
        self.timeresolution = self.simduration/self.realduration
        self.all_ID = []
        self.all_initpos = []
        self.boxsize = boxsize # has to be of the form [x,y]
        
        assert isinstance(self.boxsize, list) and len(self.boxsize) == 2, "'boxsize' must be a list of 2 elements." 
        # returns error if boxsize is not of the right format, i.e. [x,y]
        
    def Identification(self):
        for agent in range(self.popsize):
            self.all_ID = np.append(self.all_ID,agent)  
            
    def Positions(self):
        
        for ID in self.all_ID:
           self.agent_initx = rd.choice(range(self.boxsize[0]))
           self.agent_inity = rd.choice(range(self.boxsize[1]))
           
           self.all_initpos.append([self.agent_initx, self.agent_inity])
    
    def Timeline(self):
            
        for iteration in range(self.simduration):
            self.timecount += self.timeresolution
            self.timetracking = np.append(self.timetracking, self.timecount)

###----------BAT_JAMMING_00----------###
### Class implementation for agents moving:
### - within defined boundaries;
### - in a straight direction;
### - while calling at regular IPI;
### - only timing of sound taken into account;
### - without echoes.
### Initiating Inputs: ID, movdirection, iterdist, IPI
### ID: positive integer. ID of the agent.
### movdirection: integer within [-pi;pi]. Angle (in RADIANS) of the agent's direction along the plane. 
### ### mov stands for movement throughout the code.
### iterdist: positive integer. Distance covered by the agent in a single iteration. 
### ### iter & dist respectively stand for iteration & distance throughout the code.
### IPI: inter-pulse interval

###----------INITIATING OUTPUTS----------###
### boxsize: list of 2 elements. Extracted from env.boxsize (env: Launcher-class object).
### simduration: positive integer. Extracted from env.simduration (env: Launcher-class object).
### initpos: list of 2 elements. Extracted from env.all_initpos (env: Launcher-class object).
### ### init & pos respectively stand for initial & position throughout the code.
### x: positive integer. Initial abscissa of the agent
### y: positive integer. Initial ordinate of the agent
### callstarttime: positive integer. Starting time within the simulation for the agent's first call
### xhistory: list of length 1 x simduration. All abscissae of the agent through time. Initially, xhistory = x
### yhistory: list of length 1 x simduration. All ordinates of the agent through time. Initially, yhistory = y
### callshistory: array of dimensions 1 x number of calls. All calls' starting times of the agent through time.

#----------MOVEMENT----------#
# Function changing agent's position over time. Straight line with a direction set by movdirection.
# Inputs: x, y, iterdist, movdirection taken from __init___ method

#---------OUTPUTS----------#
# newx: new coordinate of the agent on the x-axis at iteration i+1
# newy: new coordinate of the agent on the y-axis at iteration i+1
# xhistory: list. updated record of agent's all past positions on the x-axis
# yhistory: list. updated record of agent's all past positions on the y-axis
#----------end of documentation for MOVEMENT----------#

#----------CALLING----------#
# Inputs: simduration, callstarttime, IPI taken from __init__ method.

#----------OUTPUTS----------#
# callshistory: array. Updated records of calls. 
# For each iteration, the output is 1 if the agent calls, and 0 otherwise.

#----------BOUNDARIES----------#
# Function keeping agent within defined boundaries; 
# when agent goes out of the boundaries, it is set back to the beginning.
# Inputs: coord, coordbound
# coord: coordinate (x or y) to be checked for i.e. position of the agent
# coordbound: boundaries for agent flow

#----------OUTPUTS----------#
# coord or Repositionned coord of the agent
#----------end of documentation for BOUNDARIES----------#

#----------SHAPEOFBEAM----------#
# IN CONSTRUCTION. Supposed to create agent's calling sonar
# Inputs: angle, shape, depth
# angle: call angle indicating width range
# shape: shape of call beam. Triangle, Hemicircle or else.
# depth: length / how far the beam can reach in frontal direction

#----------OUTPUTS----------#
# ...
#----------end of docmentation for SHAPEOFBEAM----------#
###----------end of documentation for MOVESNGROOVES----------###

class Bat_Jamming_00:
    
    def __init__(self, ID, movdirection, iterdist, IPI):
        self.ID = ID
        self.boxsize = env.boxsize
        self.simduration = env.simduration
        self.movdirection = movdirection
        self.iterdist = iterdist
        self.initpos = env.all_initpos[self.ID]
        self.x = self.initpos[0]
        self.y = self.initpos[1]        
        self.IPI = IPI
        self.callstarttime = 0 # I set it to zero at the moment, for lack of a better idea
    
        self.xhistory = [self.x]
        self.yhistory = [self.y]
        self.callshistory = np.empty([1,self.simduration], dtype = int)
        
        assert self.movdirection <= m.pi and self.movdirection >= -(m.pi), "'movdirection' must be in radians & comprised between -pi & pi."
    
    def Movement(self):
        
        for iteration in range(self.simduration):
            self.newx = self.x + self.iterdist * m.cos(self.movdirection)
            self.newy = self.y + self.iterdist * m.sin(self.movdirection)
        
            self.x = self.Boundaries(self.newx, self.boxsize[0])
            self.y = self.Boundaries(self.newy, self.boxsize[1])
        
            self.xhistory = np.append(self.xhistory, self.x)
            self.yhistory = np.append(self.yhistory, self.y)
        
    def Calling(self):
        
        for iteration in range(self.simduration):
            self.calltest = float(iteration-self.callstarttime)/float(self.IPI)
            
            if self.calltest%1 == 0:
                self.callshistory[0,iteration] = 1
            else:
                self.callshistory[0,iteration] = 0
            
            
    #def ShapeofBeam(self, angle, shape, depth): # function in construction
        #self.halfangle = angle/2
        #self.shape = shape
        #self.altitude = depth
        
        #if self.shape == "triangle":
            #self.origin = self.currpos
            #self.sidevertice = self.altitude * m.cos(self.halfangle)
            #self.range = self.sidevertice * m.cos(m.pi/2 - self.halfangle)
            #self.uppernode = [self.currpos + self.altitude, self.currpos + self.range]
            #self.lowernode = [self.currpos + self.altitude, self.currpos - self.range]
            #self.nodes = [self.currpos, self.uppernode, self.lowernode]
        #else:
            #raise AssertionError("Error: Beam shape unknown. The function ShapeofBeam covers the following shapes: 'triangle'")
            
    def Boundaries(self, coord, coordbound):
        
        if coord in range(coordbound): 
            return coord
        else:
            return range(coordbound)[0]
        print "Bat within boundaries"

### Run a multi-agents simulation and plot agents' movements.

POPSIZE = 3
BOXSIZE = [200,200]
SIMDURATION = 100
REALDURATION = 100

env = Launcher(POPSIZE, BOXSIZE, SIMDURATION, REALDURATION)
env.Identification()
env.Positions()
env.Timeline()

MOVDIRECTION = 0
ITERDIST = 10
INTER_PULSE_INTERVAL = 3

all_bats = {}
all_x = []
all_y = []

for ID in env.all_ID:
    all_bats[ID] = Bat_Jamming_00(int(ID),MOVDIRECTION,ITERDIST,INTER_PULSE_INTERVAL)
    all_bats[ID].Movement()
    plt.plot(all_bats[ID].xhistory, all_bats[ID].yhistory, marker = '^')
