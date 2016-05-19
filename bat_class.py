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
### all_ID: empty list of dimensions 1x(n=popsize elements). 
### all_initpos: empty array of n=popsize lists of 2 elements each. 
### ### init & pos respectively stand for initial and position throughout the code
### timecount: positive integer. Timing of simulation iteration.
### timeclock: array of dimension 1 x simduration. 
### timeresolution: positive float. Resolution of the time simulated.

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
# timeclock: array of dimensions 1xsimduration. Updated timeclock with "Real times" corresponding to each simulation iteration
#----------end of documentation for TIMELINE----------#
###----------end of documentation for LAUNCHER----------#  

class Launcher:
    
    def __init__(self, popsize, boxsize, simduration, realduration):
        self.popsize = popsize
        self.simduration = simduration
        self.timecount = 0 
            # set the time at 0.
        self.timeclock = np.empty(self.simduration, dtype = float)
        self.timeclock[0] = self.timecount 
            # record the first time, i.e. 0.
        self.realduration = realduration
        self.timeresolution = float(self.simduration/self.realduration) 
            # the time resolution corresponds to the relative time between simulation time and "actual" time, "percieved" by the agent.
        self.all_ID = np.empty(self.popsize, dtype = int)
        self.all_initpos = np.empty([self.popsize,2], dtype = float)
            # create an empty array of intial positions for the whole population.
            # ultimately, I think all of our array/list variables should be set this way, 
            # as it is less time-consuming computationally-wise.
        self.boxsize = boxsize 
            # has to be of the form [x,y], see below:
        
        assert isinstance(self.boxsize, list) and len(self.boxsize) == 2, "'boxsize' must be a list of 2 elements." 
            # returns error if boxsize is not of the right format, i.e. [x,y]
        
    def Identification(self):
        
        for agent in range(self.popsize):
            self.all_ID[agent] = agent
                # fills-in a list of all agents ID from 0 to n = (population size - 1)
            
    def Positions(self):
        
        for ID in self.all_ID:
           self.agent_initx = rd.choice(range(self.boxsize[0]))
               # gives each agent a random coordinate on the x-axis, within the environment boundaries
           self.agent_inity = rd.choice(range(self.boxsize[1]))
               # gives each agent a random coordinate on the y-axis, within the environment boundaries
           
           self.all_initpos[int(ID),0] = self.agent_initx
               # stores the x-coordinate in an array
           self.all_initpos[int(ID),1] = self.agent_inity
               # stores the y-coordinate in an array
    
    def Timeline(self):
            
        for iteration in range(self.simduration):
            self.timecount += self.timeresolution
                # increment timecount with timeresolution. 
                # timecount corresponds to the "actual" time, as opposed to simulated time
            self.timeclock[iteration] = self.timecount
                # store the new, incremented timecount in timeclock.

###----------BAT_JAMMING_00----------###
### Class implementation for agents moving:
### - within defined boundaries;
### - in a straight direction;
### - while calling at regular IPI;
### - only timing of sound taken into account;
### - without echoes.
### Initiating Inputs: ID, movdirection, stepsize, IPI
### ID: positive integer. ID of the agent.
### movdirection: integer within [-pi;pi]. Angle (in RADIANS) of the agent's direction along the plane. 
### ### mov stands for movement throughout the code.
### stepsize: positive integer. Distance covered by the agent in a single iteration. 
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
# Inputs: x, y, stepsize, movdirection taken from __init___ method

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
    
    def __init__(self, ID, movdirection, stepsize, IPI):
        self.ID = ID
        self.boxsize = env.boxsize
            # takes boxsize from env, an object of the class Launcher
        self.simduration = env.simduration
            # takes simduration from env, an object of the class Launcher
        self.movdirection = movdirection
        self.stepsize = stepsize       
        self.IPI = IPI
        
        assert self.movdirection <= m.pi and self.movdirection >= -(m.pi), "'movdirection' must be in radians & comprised between -pi & pi."
            # returns an error message if movdirection is not within [-pi;pi].
    
    def Movement(self):
        self.newx = self.x + self.stepsize * m.cos(self.movdirection)
            # calculates the new x coordinate according to:
                # the distance travelled over 1 time step.
                # the direction of the movement
        self.newy = self.y + self.stepsize * m.sin(self.movdirection)
             # calculates the new y coordinate according to:
                # the distance travelled over 1 time step.
                # the direction of the movement
        
        self.x = self.Boundaries(self.newx, self.boxsize[0])
            # verify that the new x coordinate is within boundaries
        self.y = self.Boundaries(self.newy, self.boxsize[1])
            # verify that the new y coordinate is within boundaries
        
        self.xhistory = np.append(self.xhistory, self.x)
            # store the newly calculated x in xhistory
        self.yhistory = np.append(self.yhistory, self.y)
            # store the newly calculated y in yhistory
        
    def Calling(self):
        self.calltest = float(self.timestep-self.callstarttime)/float(self.IPI)
            # calculates the theoretical number of calls since the 1st call:
            # (time since the 1st call = current time - time of the 1st call)/IPI
            
        if self.calltest%1 == 0:
            self.callshistory = np.append(self.callshistory, 1)
                # adds a new call (accounted for with a 1) if the theoretical number of 
                # calls since the 1st call is a natural number
        else:
            self.callshistory = np.append(self.callshistory, 0)
                # adds a no new call (accounted for with a 0) if the theoretical number
                # of calls since the 1st call is a float            
            
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
                # if the coordinate is within boundaries, it stays the same.
        else:
            return range(coordbound)[0]
                # if the coordinate isn't within boundaries, it is reset to the beginning.
        print "Bat within boundaries"

### Run a multi-agents simulation and plot agents' movements.

POPSIZE = 3
    # size of the population of bats you want to simulate
BOXSIZE = [200,200]
    # space within which the bats are moving
SIMDURATION = 100
    # duration of the simulation you want to run
REALDURATION = 100 
    # duration of the simulation in "real" time
    # e.g. you can run the simulation for say,
    # 100 iterations, corresponding to a total of 1000 milliseconds.

env = Launcher(POPSIZE, BOXSIZE, SIMDURATION, REALDURATION)
    # sets the simulation environment with Launcher, according to the parameters given 
    # above, and stores it in an object called env
env.Identification()
env.Positions()
env.Timeline()
    # run the functions contained in the Launcher class

MOVDIRECTION = 0
    # direction of the movement of the bat. 
STEPSIZE = 10
    # distance travelled by the bat over each time step / iteration
INTER_PULSE_INTERVAL = 3
    # Inter-pulse interval, i.e. time interval between each call 
    # For now these 3 parameters are the same for all agents, 
    # but this can be changed if we want to

all_bats = {}
    # creates an empty dictionary for every bat instance to be stored
all_x = {key:{0:0} for key in list(env.all_ID)}
all_y = {key:{0:0} for key in list(env.all_ID)}

for ID in env.all_ID:
    all_bats[ID] = Bat_Jamming_00(int(ID), MOVDIRECTION,STEPSIZE,INTER_PULSE_INTERVAL)
        # stores all instances of the class Bat_Jamming_00 within the bat population
    all_bats[ID].x = env.all_initpos[int(ID)][0]
        # initial x coordinate for each instance, taken from env
    all_bats[ID].xhistory = [all_bats[ID].x]
        # stores it in xhistory
    all_bats[ID].y = env.all_initpos[int(ID)][1]
        # initial y coordinate for each instance, taken from env
    all_bats[ID].yhistory = [all_bats[ID].y]
        # stores it in yhistory
    all_bats[ID].callstarttime = 0 
        # sets the starting time for the initial call 
        # set to zero at the moment, for lack of a better idea
    all_bats[ID].callshistory = []
        # creates an empty list to store calls records
    all_x[ID][0] = all_bats[ID].xhistory[0]
    all_y[ID][0] = all_bats[ID].yhistory[0]


for timestep in env.timeclock:
    for ID in env.all_ID:
        all_bats[ID].timestep = timestep
            # indicates the current time step for each instance 
        all_bats[ID].Movement()
            # makes the instance move
        
        if ID not in all_x:
            all_x.add(ID)
            
            if timestep not in all_x:
                all_x.update({ID:{timestep:all_bats[ID].xhistory[timestep]}})
            else:
                continue
        
        else:
            continue
        
        all_y.update({ID:{timestep:all_bats[ID].yhistory[timestep]}})
        all_bats[ID].Calling()
            # makes the instance call
        plt.plot(all_bats[ID].xhistory, all_bats[ID].yhistory, marker = '^')
            # plot all instances movements over time

timestep = 2
ID = 1
callsources = {}
max_timestore = 3

def Min_circle(center_x, center_y, radius, x, y):
    # nb: x and y have to be np.arrays for this function to work
    dist = (x - center_x) ** 2 + (y - center_y) ** 2
    return dist > radius

def Max_circle(center_x, center_y, radius, x, y):
    # nb: x and y have to be np.arrays for this function to work
    dist = (x - center_x) ** 2 + (y - center_y) ** 2
    return dist < radius
    
def TOA(ID, timestep):
    
    if all_bats[ID].callshistory[timestep] == 1:
        callsources.update({ID:{timestep:{'xsource': all_bats[ID].xhistory[timestep], 'ysource': all_bats[ID].yhistory[timestep], 'propdist':0}}})
        # if this bat called at this timestep, store the position of the bat at the time 
        # of calling into a dictionary of the form {bat:{time of calling:[x,y,propdist]}}
        # I have to check whether .update is actually the appropriate way to do it        
        
    else:
        continue
    
    if ID in callsources.keys():
    # if the agent has ever called before: 
        
        for calltime in callsources[ID].keys():
        # for each time step at which the agent previously called:
            timestore = timestep - calltime
            # the storing time of the corresponding call into the dictionary 
            # is calculated. Remember that iteration time and "real" time can 
            # be different --> choose max_timestore appropriately
            
            if timestore > max_timestore:
            # if the storing time is longer than a predefined time: 
                callsources.pop(callsources[ID].keys()[calltime], None)
                # erase it from the memory / dictionary
            
            else:
                continue
            
            speed_sound = 340.29 # speed of sound at sea level in m/s
            propdist = float(speed_sound * timestore / env.timeresolution)
            # calculate, for a certain call, its propagation distance at timestep 
            # according to the time when the call was emitted and the speed of sound
            callsources[ID][timestep]['propdist'] = propdist
            # update the propagation distance in callsource
        
    else:
        continue
    
            
            