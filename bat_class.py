# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 14:31:31 2016

@author: Claire
"""

import random as rd
import math  as m
import numpy as np
import matplotlib.pyplot as plt

Duration_of_Simulation = 100
Size_of_Box = [200,200]
Size_of_Population = 100

###----------LAUNCHER----------###
### Class implementation to launch the environment for bats simulations
### Initiating inputs: simduration, boxsize, x, y
### simduration: positive integer. Duration of simulation 
### ### sim stands for simulation throughout code
### boxsize: list of 2 elemnts. Cartesian space within which the agents can move
### all_x: list of n elements, where n is the number of agents. Initial abscissa(e) of agent(s)
### all_y: list of n elements, where n is the number of agents. Initial ordinate(s) of agent(s)

###----------INITIATING OUTPUTS----------###
### all_initpos: list of 2 lists of n elements each. Initial position(s) of agent(s)
### ### init & pos respectively stand for initial and position throughout the code
### timecount: positive integer. Timing of the simulation

#----------TIMELINE----------#
# Function drawing the timeline of the simulation
# Inputs: taken from __init__ method

#----------OUTPUTS----------#
# timecount: updated timecount. At the end of the simulation, timecount = simduration
# timetracking: array of dimensions 1 x simduration. Timeline of the simulation, with size 1 timesteps.
#----------end of documentation for TIMELINE----------#
###----------end of documentation for LAUNCHER----------#  

class Launcher:
    
    def __init__(self, popsize, simduration, boxsize):
        self.popsize = popsize
        self.simduration = simduration
        self.timecount = 0
        self.all_ID = []
        self.all_initpos = []
        self.boxsize = boxsize # has to be of the form [x,y]
        
        assert isinstance(self.boxsize, list) and len(self.boxsize) == 2, "'boxsize' must be a list of 2 elements." 
        # returns error if boxsize is not of the right format, i.e. [x,y]
        
    def Identification(self):
        for i in range(self.popsize):
            self.all_ID = np.append(self.all_ID,i)    
    
    def Timeline(self):
        self.timetracking = [self.timecount]
    
        for iteration in range(self.simduration):
            self.timecount += 1
            self.timetracking = np.append(self.timetracking, self.timecount)
    
    def Positions(self):
        
        for ID in self.all_ID:
           self.agent_initx = rd.choice(range(self.boxsize[0]))
           self.agent_inity = rd.choice(range(self.boxsize[1]))
           self.agent_initpos = [self.agent_initx, self.agent_inity]
           
           self.all_initpos.append(self.agent_initpos)
            
env = Launcher(Size_of_Population, Duration_of_Simulation, Size_of_Box)
env.Identification()
env.Positions()
env.Timeline()


###----------BAT_JAMMING----------###
### Class implementation for agent moving within defined boundaries, in a straight direction.
### Initiating Inputs: movdirection, iterdist, simduration
### movdirection: integer within [-pi;pi]. Angle (in RADIANS) of the agent's direction along the plane. 
### ### mov stands for movement throughout the code.
### iterdist: positive integer. Distance covered by the agent in a single iteration. 
### ### iter and dist respectively stand for iteration & distance throughout the code.
### simduration: positive integer. Duration of simulation.
### ### sim stands for simulation throughout the code.
### IPI: inter-pulse interval
### callduration: timing of call
### boxsize: space within which the agent can move 

###----------INITIATING OUTPUTS----------###
### x: positive integer. Initial abscissa of the agent, extracted from all_initpos (Launcher class)
### y: positive integer. Initial ordinate of the agent, extracted from all_initpos (Launcher class)
### initpos: vector of length 2. Position of the agent (comprises x and y)
### ### init & pos respectively stand for initial & position throughout the code.
### callstarttime: positive integer. Starting time within the simulation for the agent's first call
### callendtime: positive integer. Ending time of the call.
### positionshistory: list of length 2*simduration. All positions of the agent through time
### callstiminghistory: list of length 2*number of calls. All calls' starting and ending time of the agent through time.

#----------BOUNDARIES----------#
# Function keeping agent within defined boundaries; 
# when agent goes out of the boundaries, it is set back to the beginning.
# Inputs: coord, coordbound
# coord: coordinate (x or y) to be checked for i.e. position of the agent
# coordbound: boundaries for agent flow

#----------OUTPUTS----------#
# coord or Repositionned coord of the agent
#----------end of documentation for BOUNDARIES----------#

#----------MOVEMENT----------#
# Function changing agent's position over time. Straight line when turnangle = 0.
# Inputs: taken from __init___ OUTPUTS

#---------OUTPUTS----------#
# newx: new coordinate of the agent on the x-axis
# newy: new coordinate of the agent on the y-axis
# newpos: new coordinates of the agents at iteration i+1 (x, y)
# positionhistory: list. updated record of agent's all past positions
#----------end of documentation for MOVEMENT----------#

#----------PLOTTING----------#
# Function compiling graphic visualisation of agent's movement over time.
# Inputs: taken from MOVEMENT function OUTPUTS

#----------OUTPUTS----------#
# xAxis: vector of all x
# yAxis: vector of all y
#----------end of documentation for PLOTTING----------#

#----------CALLING----------#
# Inputs: taken from __init__ method.

#----------OUTPUTS----------#
# callhistory

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

Direction_of_Movement = 0
Dist_Covered_Per_Iter = 10
Inter_Pulse_Interval = 3

class Bat_Jamming:
    
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
        self.callstarttime = 0 # I set it at zero at the moment, for lack of a better idea
    
        self.positionshistory = [self.x, self.y]
        self.callshistory = np.empty([1,self.simduration], dtype = int)
        
        assert self.movdirection <= m.pi and self.movdirection >= -(m.pi), "'movdirection' must be in radians & comprised between -pi & pi."
    
    def Movement(self):
        
        for iteration in range(self.simduration):
            self.newx = self.x + self.iterdist * m.cos(self.movdirection)
            self.newy = self.y + self.iterdist * m.sin(self.movdirection)
        
            self.x = self.Boundaries(self.newx, self.boxsize[0])
            self.y = self.Boundaries(self.newy, self.boxsize[1])
        
            self.currpos = [self.x, self.y]
        
            self.positionshistory = np.append(self.positionshistory, [self.x, self.y])
        
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
        
    def Plotting(self):
        self.xAxis = [] 
        self.yAxis = []
        
        for x, y in zip(self.positionshistory[0:2*self.simduration:2],self.positionshistory[1:2*self.simduration:2]):
            self.xAxis = np.append(self.xAxis, x)
            self.yAxis = np.append(self.yAxis, y)
        
        plt.plot(self.xAxis, self.yAxis, 'ro')