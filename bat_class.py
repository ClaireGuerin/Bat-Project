# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 14:31:31 2016

@author: Claire
"""

import random as rd
import math  as m
import numpy as np
import matplotlib.pyplot as plt
import pylab

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

###----------end of documentation for BAT_JAMMING_00----------###

class Bat_Jamming_00:
    
    def __init__(self, ID, movdirection, stepsize, IPI, max_timestore, ring_dev):
        self.ID = int(ID)
        self.boxsize = env.boxsize
            # takes boxsize from env, an object of the class Launcher
        self.simduration = int(env.simduration)
            # takes simduration from env, an object of the class Launcher
        self.movdirection = int(movdirection)
        self.stepsize = int(stepsize)       
        self.IPI = int(IPI)
        self.max_timestore = max_timestore
        self.ring_dev = int(ring_dev)
        self.callsource = {}
        self.hearhistory = []
        
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
                # adds a "no new call" (accounted for with 0) if the theoretical number
                # of calls since the 1st call is a float            
            
    def TOA(self, timestep):
    
        if all_bats[self.ID].callshistory[timestep] == 1:
            Dict_update(self.callsource, {self.ID:{timestep:{'xsource': all_bats[ID].xhistory[timestep], 'ysource': all_bats[ID].yhistory[timestep], 'propdist':0}}})
                # if this bat called at this timestep, store the position of the bat at the time 
                # of calling into a dictionary of the form {bat:{time of calling:[x,y,propdist]}}
                # I have to check whether D_update is actually the appropriate way to do it        
    
        if self.ID in self.callsource.keys():
            # if the agent has ever called before: 
        
            for calltime in self.callsource[self.ID].keys():
                # for each time step at which the agent previously called:
                self.Data_storage(self.callsource, timestep, calltime)
                    # erase calls that have been emitted too long ago to be heared anymore
                self.Propagation(self.callsource, calltime)
                    # update the propagation distance in callsource
            
                for identity in env.all_ID:
                    # for each agent in the simulation:
                
                    if identity in self.callsource.keys():
                        # if the agent has ever called before: 
                    
                        for calltime in self.callsource[self.ID].keys():
                            # for each time step at which the agent previously called:
                            self.Hearing_test(identity, calltime, self.hearhistory, timestep)
                                # identify and record calls that can be heard.
    
    def Boundaries(self, coord, coordbound):
        
        if coord in range(coordbound): 
            return coord
                # if the coordinate is within boundaries, it stays the same.
        else:
            return range(coordbound)[0]
                # if the coordinate isn't within boundaries, it is reset to the beginning.
        print "Bat within boundaries"
        
    def Data_storage(self, dict1, tmstp, tcall):
        self.timestore = tmstp - tcall
            # storing time of the corresponding call into the dictionary 
            # Remember that iteration time and "real" time can be different 
            # --> choose max_timestore appropriately
            
        if self.timestore > self.max_timestore:
            # if the storing time is longer than a predefined time: 
            dict1[self.ID].pop(tcall, None)
                # erase it from the memory / dictionary
        
        return dict1

    def Min_circle(self, center_x, center_y, radius, x, y):
        dist = (x - center_x) ** 2 + (y - center_y) ** 2
        return dist > radius

    def Max_circle(self, center_x, center_y, radius, x, y):
        dist = (x - center_x) ** 2 + (y - center_y) ** 2
        return dist < radius
        
    def Propagation(self, dict1, tmstp):
        SPEED_SOUND = 340.29 # speed of sound at sea level in m/s
        self.propdist = float(SPEED_SOUND * self.timestore / env.timeresolution)
            # calculate, for a certain call, its propagation distance at timestep 
            # according to the time when the call was emitted and the speed of sound
        dict1[self.ID][tmstp]['propdist'] = self.propdist
            # update the propagation distance in dict1
        
        return dict1
    
    def Hearing_test(self, ag_id, tcall, nparray, tmstp):
        self.mintest = self.Min_circle(self.callsource[ag_id][tcall]['xsource'],
                                       self.callsource[ag_id][tcall]['ysource'],
                                       self.callsource[ag_id][tcall]['propdist'] - self.ring_dev,
                                       all_bats[int(ID)].xhistory[tmstp],
                                       all_bats[int(ID)].yhistory[tmstp])
                                         
        self.maxtest = self.Max_circle(self.callsource[ag_id][tcall]['xsource'], 
                                       self.callsource[ag_id][tcall]['ysource'], 
                                       self.callsource[ag_id][tcall]['propdist'] + self.ring_dev,
                                       all_bats[int(ID)].xhistory[tmstp],
                                       all_bats[int(ID)].yhistory[tmstp])
                
        if self.mintest and self.maxtest:
            # test if the agent 'ID' can hear any call from agent 'identity'
            # including own calls
            nparray = np.append(nparray,[self.ID, tmstp, ag_id, tcall])
                # if so, store the IDs of the hearing bat and the calling bat,
                # as well as the time at which the bat called and has been heard
        
        return nparray

###----------HEARING_00----------###
### Class implementation for agents hearing:
### - calls from the others and themselves;
### - considering sound propagates in a ring shape, towards all directions
### Initiating inputs: ID, timestep, max_timestore
### ID: integer. Focal bat/agent identification number for which sound hearing is recorded
### timestep: integer. Simulation time step at which the bat is hearing
### max_timestore: integer. Maximum time for storing the data, in simulation steps.
###     Corresponds to the time after which the intensity of the call is considered too low to be heard anymore.
###     Ideally, should implement hearing threshold instead e.g. 20 (pe)SPL

###----------INITIATING OUTPUTS----------###
### callsources: dictionary. Keeps track of the coordinates of all calls emitted by 
###     every agent, at each timestep, as well as the updated propagation distance from the 
###     source at the current time
### hearhistory: numpy array of dimensions n * 4. Reports when an agent heard a sound.
###     Data provided is, in order: 
###     - ID of the hearing agent
###     - time of hearing
###     - ID of the calling agent = source identification
###     - time at which the source called

#----------TOA----------#
# Function recording the calls heard by a specific agent, at a soecific time.
# Inputs: ID, timestep, max_timestore, callsources, hearhistory taken from __init__ method

#----------OUTPUTS----------#
# callsources: dictionary. Updated record of all calls sources
# callshearing: numpy array. Updated record of all heard calls
#----------End of documentation for TOA----------#

#----------Dict_update----------#
# Function for updating a dictionary without over-writing the keys already stored
# in the dictionary.
# Inputs: dict1, dict2
# dict1: original dictionary to be updated.
# dict2: dictionary to add to dict1. NB: has to be of the same format as dict1.

#----------OUTPUTS----------#
# dict1: updated dictionary (new key is added to the already existing keys in the 
#   dictionary, instead of replacing them) 
#----------End of documentation for DICT_UPDATE----------#

#----------DATA_STORAGE----------#
# Function applying a time limit for storing keys in dictionary 
# It is meant for callsources in TOA function, each sources must be deleted when the 
# call is supposedly not likely to be heard anymore for it was emitted too long ago to 
# show a significant intensity of call
# Inputs taken from __init__ method: ID, timestep, max_timestore
# Inputs: dict1, tcall
# dict1: dictionary to be scrutinized for "out-of-date" calls. Here, dict1 will be callsources
# tcall: integer. Time at which the call was emitted

#----------OUPUTS----------#
# timestore: integer. Time since when the call information has been stored in the dictionary
#   (callsources) in timesteps of simulation
# dict1: dictionary. Updated such as tcalls > maximum time of storage have been deleted
#----------End of documentation for DATA_STORAGE----------#  

#----------MIN_CIRCLE----------#
# Function defining the inside concentric circle of the ring of propagation of a call
# and checking if a specific point is outside the circle. Returns a boolean.
# NB: x and y have to be int or np.arrays for this function to work
# Inputs: center_x, center_y, radius, x, y 
# center_x: integer. Abscissa of the center of the circle (i.e. source)
# center_y: integer. Ordinate of the center of the circle (i.e. source)
# radius: integer. radius of the circle
# x: integer. Abscissa of the point to be checked for
# y: integer. Ordinate of the point to be checked for

#----------OUTPUTS----------#
# dist: integer. Distance of the specific point to the center of the circle.
# If dist is more than radius, the function returns TRUE.
#----------End of documentation for MIN_CIRCLE----------#

#----------MAX_CIRCLE----------#
# Function defining the outside concentric circle of the ring of propagation of a call
# and checking if a specific point is inside the circle. Returns a boolean.
# NB: x and y have to be int or np.arrays for this function to work
# Inputs: center_x, center_y, radius, x, y 
# center_x: integer. Abscissa of the center of the circle (i.e. source)
# center_y: integer. Ordinate of the center of the circle (i.e. source)
# radius: integer. radius of the circle
# x: integer. Abscissa of the point to be checked for
# y: integer. Ordinate of the point to be checked for

#----------OUTPUTS----------#
# dist: integer. Distance of the specific point to the center of the circle.
# If dist is less than radius, the function returns TRUE.
#----------End of documentation for MAX_CIRCLE----------#

#----------PROPAGATION----------#
# Function for calculation of the distance of the sound from its source as it propagates
# Inputs taken from __init__ method: ID, timestep
# Inputs: dict1
# Other inputs: timestore taken from data_storage method
# dict1: dictionary (here, callsource) from which the propagation distance is extracted

#----------OUTPUTS----------#
# propdist: integer. Distance of propagation of the sound from its source
# dict1: dictionary. callsource dictionary with updated propagation distances
#----------End of documentation for PROPAGATION----------#

#----------HEARING_TEST----------#
# Function that checks, for a given call, whether the focal bat/agent is able to hear it
# The agent can hear a call if and only if both conditions Min_circle and Max_circle are
# met. In other words, if the focal bat is located in between the two circles determining
# the ring of sound propagated from a call source at a specific time.
# Inputs taken from __init__ method: callsources, ring_dev, timestep, ID
# Inputs: ag_id, tcall, nparray
# ag_id: integer. Identification number of the focal agent
# tcall: integer. time at which a given call was emitted
# nparray: numpy array. Meant to be hearhistory

#----------OUTPUTS----------#
# mintest: boolean. Is the bat out of the inside circle of the ring of sound?
# maxtest: boolean. Is the bat in the outside circle of the ring of sound?
# nparray: numpy array. Updated (hearhistory) with the call and hearing times, if the 
# call is indeed heard by the focal agent/bat
#----------End of documentation for HEARING_TEST----------#

###----------End of documentation for HEARING_00----------###  
       
class Hearing_00:
    
    def __init__(self, ID, timestep, max_timestore, ring_dev):
        self.max_timestore = int(max_timestore)
        self.ID = int(ID)
        self.timestep = int(timestep)
        self.ring_dev = int(ring_dev)
        self.callsource = {}
        self.hearhistory = []
    
    def TOA(self):
    
        if all_bats[self.ID].callshistory[self.timestep] == 1:
            self.callsource = {self.ID:{self.timestep:{'xsource': all_bats[ID].xhistory[self.timestep], 'ysource': all_bats[ID].yhistory[self.timestep], 'propdist':0}}}
                # if this bat called at this timestep, store the position of the bat at the time 
                # of calling into a dictionary of the form {bat:{time of calling:[x,y,propdist]}}
                # I have to check whether D_update is actually the appropriate way to do it        
    
        if self.ID in self.callsource.keys():
            # if the agent has ever called before: 
        
            for calltime in self.callsource[self.ID].keys():
                # for each time step at which the agent previously called:
                self.Data_storage(self.callsource, calltime)
                    # erase calls that have been emitted too long ago to be heared anymore
                self.Propagation(self.callsource)
                    # update the propagation distance in callsource
            
                for identity in env.all_ID:
                    # for each agent in the simulation:
                
                    if identity in self.callsource.keys():
                        # if the agent has ever called before: 
                    
                        for calltime in self.callsource[self.ID].keys():
                            # for each time step at which the agent previously called:
                            self.Hearing_test(identity, calltime, self.hearhistory)
                                # identify and record calls that can be heard.

    def Data_storage(self, dict1, tcall):
        self.timestore = self.timestep - tcall
            # storing time of the corresponding call into the dictionary 
            # Remember that iteration time and "real" time can be different 
            # --> choose max_timestore appropriately
            
        if self.timestore > self.max_timestore:
            # if the storing time is longer than a predefined time: 
            dict1.pop(dict1[self.ID].keys()[tcall], None)
                # erase it from the memory / dictionary
        
        return dict1

    def Min_circle(self, center_x, center_y, radius, x, y):
        dist = (x - center_x) ** 2 + (y - center_y) ** 2
        return dist > radius

    def Max_circle(self, center_x, center_y, radius, x, y):
        dist = (x - center_x) ** 2 + (y - center_y) ** 2
        return dist < radius
        
    def Propagation(self, dict1):
        SPEED_SOUND = 340.29 # speed of sound at sea level in m/s
        self.propdist = float(SPEED_SOUND * self.timestore / env.timeresolution)
            # calculate, for a certain call, its propagation distance at timestep 
            # according to the time when the call was emitted and the speed of sound
        dict1[self.ID][self.timestep]['propdist'] = self.propdist
            # update the propagation distance in dict1
        
        return dict1
    
    def Hearing_test(self, ag_id, tcall, nparray):
        self.mintest = self.Min_circle(self.all_sources[ag_id][tcall]['xsource'],
                                       self.all_sources[ag_id][tcall]['ysource'],
                                       self.all_sources[ag_id][tcall]['propdist'] - self.ring_dev,
                                       all_bats[int(ID)].xhistory[self.timestep],
                                       all_bats[int(ID)].yhistory[self.timestep])
                                         
        self.maxtest = self.Max_circle(self.all_sources[ag_id][tcall]['xsource'], 
                                       self.all_sources[ag_id][tcall]['ysource'], 
                                       self.all_sources[ag_id][tcall]['propdist'] + self.ring_dev,
                                       all_bats[int(ID)].xhistory[self.timestep],
                                       all_bats[int(ID)].yhistory[self.timestep])
                
        if self.mintest and self.maxtest:
            # test if the agent 'ID' can hear any call from agent 'identity'
            # including own calls
            nparray = np.append(nparray,[self.ID, self.timestep, ag_id, tcall])
                # if so, store the IDs of the hearing bat and the calling bat,
                # as well as the time at which the bat called and has been heard
        
        return nparray



def Dict_update(dict1, dict2):
    
    for key in dict2:
        
        if key in dict1:
            dict1[key].update(dict2[key])
        else:
            dict1[key] = dict2[key]
            
    return dict1
        
### Run a multi-agents simulation and plot agents' movements.
###----------PARAMETERS----------###
### to be determined / chosen beforehand 
### POPSIZE: size of the population of bats you want to simulate
### BOXSIZE: space within which the bats are moving
### SIMDURATION: duration of the simulation you want to run
### REALDURATION: duration of the simulation in "real" time
###     e.g. you can run the simulation for say,
###     100 iterations, corresponding to a total of 1000 milliseconds.
### MOVDIRECTION: direction of the movement of the bat.
### STEPSIZE: distance travelled by the bat over each time step / iteration
### INTER_PULSE_INTERVAL: Inter-pulse interval, i.e. time interval between each call 
###     For now these 3 parameters are the same for all agents, but can be set as different
### MAX_TIMESTORE: maximum time (in simulation time steps) for a call source to be taken into
###     account for the hearing part of the simulation

POPSIZE = 3
BOXSIZE = [200,200]
SIMDURATION = 100
REALDURATION = 100
MOVDIRECTION = 0
STEPSIZE = 10
INTER_PULSE_INTERVAL = 3
MAX_TIMESTORE = 50 # ideally, should implement hearing threshold instead e.g. 20 (pe)SPL
RING_DEV = 50

# Set the simulation environment with Launcher, according to the parameters given 
# above, and stores it into an object called env
# Run the functions contained in the Launcher class
env = Launcher(POPSIZE, BOXSIZE, SIMDURATION, REALDURATION)
env.Identification()
env.Positions()
env.Timeline()

all_bats = {}
    # creates empty dictionaries for every bat instance to be stored
#all_x = {key:{0:0} for key in list(env.all_ID)}
#all_y = {key:{0:0} for key in list(env.all_ID)}
all_sources = {}

for ID in env.all_ID:
    all_bats[int(ID)] = Bat_Jamming_00(int(ID), MOVDIRECTION,STEPSIZE,INTER_PULSE_INTERVAL, MAX_TIMESTORE, RING_DEV)
        # stores all instances of the class Bat_Jamming_00 within the bat population
    all_bats[int(ID)].x = env.all_initpos[int(ID)][0]
        # initial x coordinate for each instance, taken from env
    all_bats[int(ID)].xhistory = [all_bats[int(ID)].x]
        # stores it in xhistory
    all_bats[int(ID)].y = env.all_initpos[int(ID)][1]
        # initial y coordinate for each instance, taken from env
    all_bats[int(ID)].yhistory = [all_bats[int(ID)].y]
        # stores it in yhistory
    all_bats[int(ID)].callstarttime = 0 
        # sets the starting time for the initial call 
        # set to zero at the moment, for lack of a better idea
    all_bats[int(ID)].callshistory = []
        # creates an empty list to store calls records
    #all_x[int(ID)][0] = all_bats[int(ID)].xhistory[0]
    #all_y[int(ID)][0] = all_bats[int(ID)].yhistory[0]


for timestep in env.timeclock:
    for ID in env.all_ID:
        all_bats[int(ID)].timestep = timestep
            # indicates the current time step for each instance 
        all_bats[int(ID)].Movement()
            # makes the instance move
        all_bats[int(ID)].Calling()
            # makes the instance call  
        all_bats[int(ID)].TOA(int(timestep-1))
        # need to find a way to have callsource "in-&-out" of the class... 
        
        Dict_update(all_sources,all_bats[int(ID)].callsource) 
         
        plt.plot(all_bats[int(ID)].xhistory, all_bats[int(ID)].yhistory, marker = '^')
            # plot all instances movements over time