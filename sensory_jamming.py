# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 14:31:31 2016

@author: Claire
"""

from __future__ import division
import random as rd
import math  as m
import numpy as np
import matplotlib.pyplot as plt

# import pylab

class Launcher:
    # Class implementation to initiate simulatory environment.

    def __init__(self, popsize, boxsize, d_t, simduration):
        self.popsize = int(popsize)
            # integer. Size of the bat population / Number of agents 
        self.d_t = float(d_t * m.pow(10,-3))
            # float. Time resolution, i.e. how much time (in ms) 
            # is expressed in 1 simulation time step.
        self.simduration = int(simduration)
            # integer. Total number of time steps to be run in the simulation. 
        self.realtime = 0 
            # time of the simulation, in seconds. Initially set at 0.
        self.timeclock = np.empty(self.simduration, dtype = float)
            # empty array for floats, dimensions 1*simduration. Will 
            # contain all times in seconds for each simulation time step.
        self.timeclock[0] = float(self.realtime) 
            # record the first time, i.e. 0.
        self.realduration = float(self.simduration * self.d_t)
            # duration of the simulation in seconds.
        self.all_ID = np.empty(self.popsize, dtype = int)
            # empty array for integers, dimensions 1*popsize. Will 
            # contain the identification numbers of all the agents to be run.
        self.all_initpos = np.empty([self.popsize,2], dtype = float)
            # empty array of intial positions for the whole population.
        self.callsources = {}
            # empty dictionary. Will contain the sources of each calls, 
            # emitted by every agent throughout the simulation.
        self.boxsize = boxsize 
            # list of 2 integers. Area (rectangle) within which the agents can 
            # move, defined by the lengths of edges in meters. 
            # Has to be of the form [x,y], which is asserted below:
        
        assert isinstance(self.boxsize, list) and len(self.boxsize) == 2, "'boxsize' must be a list of 2 elements." 
            # returns error if boxsize is not of the right format, i.e. [x,y]
        
    def Identification(self):
        
        for agent in range(self.popsize):
            # for each agent in the bat population:
            self.all_ID[int(agent)] = agent
                # fill-in a list of all agents ID 
                # from 0 to n = (population size - 1)
            
    def Positions(self):
        
        for ID in self.all_ID:
            # for each agent in the bat population:
           self.agent_initx = rd.uniform(0, self.boxsize[0])
               # give each agent a random coordinate on the x-axis, 
               # within the environment boundaries
           self.agent_inity = rd.uniform(0, self.boxsize[1])
               # give each agent a random coordinate on the y-axis, 
               # within the environment boundaries
           self.all_initpos[int(ID),0] = self.agent_initx
               # store the x-coordinate in an array
           self.all_initpos[int(ID),1] = self.agent_inity
               # store the y-coordinate in an array
    
    def Timeline(self):
            
        for iteration in range(self.simduration):
            # for each time step in the simulation:
            self.realtime += float(self.d_t)
                # increment timecount with time resolution. 
            self.timeclock[int(iteration)] = self.realtime
                # store the new, incremented time count in timeclock.

    class Bat_Jamming_00:
        # Class implementation nested into Launcher class implementation.
        # Sets an agent's movement, calling and hearing rules.
    
        def __init__(self, ID, movdirection, flightspeed, IPI, max_hear_dist):
            self.ID = int(ID)
                # integer. Identification number of the focal agent.
            self.boxsize = env.boxsize
                # takes boxsize from env (class Launcher).
            self.d_t = env.d_t
                # takes d_t from env (class Launcher).
            self.movdirection = float(movdirection)
                # float. Direction (in radians) towards which the agent flies.
                # Must be between -pi and pi (asserted later).
            self.stepsize = float(flightspeed * self.d_t) 
                # float. Distance in meters covered by the agent in 1 time step
            self.IPI = int(float(IPI * m.pow(10,-3)) / env.d_t)
                # int. Agent's inter-pulse interval, converted in time steps.
            self.speedsound = 340.29 
                # float. Speed of sound at sea level in m/s.
            self.max_timestore = float(max_hear_dist / self.speedsound)
                # float. Time maximal, in seconds, during which a sound can 
                # travel, before its intensity passes below the hearing 
                # threshold of the agent.
            self.ring_width = float(1) 
                # float. Spatial difference between start and end of call (m).
            self.hearhistory = []
                # empty array, dimensions 3*unknown. Will store the 
                # identification number of the agent who's call has been heard 
                # by the focal bat; the time at which the call was emitted; &
                # the time at which it has been heard by the focal agent.
        
            assert self.movdirection <= m.pi and self.movdirection >= -(m.pi), "'movdirection' must be in radians & comprised between -pi & pi."
                # returns an error message if movdirection is not within [-pi;pi].
    
        def Movement(self):
            self.newx = float(self.x + self.stepsize * m.cos(self.movdirection))
                # calculate the new x coordinate according to:
                # - the distance travelled over 1 time step.
                # - the direction of the movement
            self.newy = float(self.y + self.stepsize * m.sin(self.movdirection))
                # calculate the new y coordinate according to:
                # - the distance travelled over 1 time step.
                # - the direction of the movement
        
            self.x = self.Boundaries(self.newx, self.boxsize[0])
                # verify that the new x coordinate is within boundaries
            self.y = self.Boundaries(self.newy, self.boxsize[1])
                # verify that the new y coordinate is within boundaries
        
            self.xhistory = np.append(self.xhistory, self.x)
                # store the newly calculated x in xhistory
            self.yhistory = np.append(self.yhistory, self.y)
                # store the newly calculated y in yhistory
        
        def Calling(self):
            
            self.calltest = float(self.timestep-self.firstcall)/float(self.IPI)
                # calculate the theoretical number of calls since the 1st call:
                # time since 1st call = current time - time of 1st call/IPI
            
            if self.calltest%1 == 0:
                # if the theoretical number of calls since the 1st call is a 
                # natural number:
                self.callshistory[int(self.timestep)] = 1
                    # add a new call (accounted for with a 1) 
                Dict_update(env.callsources, {self.ID:{self.timestep:{'xsource': self.xhistory[self.timestep], 'ysource': self.yhistory[self.timestep], 'propdist':0}}})
                    # store the position of the bat at the time of calling into 
                    # a dictionary of the form: 
                    # {bat:{time of calling:[x,y,propdist]}}
                
            else:
                self.callshistory[int(self.timestep)] = 0
                    # add a "no new call" (accounted for with 0)


            if self.ID in env.callsources.keys():
                # if the agent has ever called before: 
        
                for calltime in env.callsources[int(self.ID)].keys():
                    # for each time step at which a call was emitted:
                    self.Data_storage(env.callsources, calltime)
                        # erase calls that are too old to be heared anymore
                    if calltime in env.callsources[int(self.ID)].keys():
                        # if the call is still in the dictionary:
                        self.Propagation(env.callsources, calltime)
                            # update its propagation distance accordingly

                    
        def Hearing(self):
            
            for identity in env.all_ID:
                # for each agent in the simulation:
                
                if identity in env.callsources.keys():
                    # if the agent has ever called before: 
                    
                    for calltime in env.callsources[int(identity)].keys():
                        # for each time step at which the agent previously called:
                        self.Hearing_test(identity, calltime, self.hearhistory)
                        # identify and record calls that can be heard by focal agent
    
        def Boundaries(self, coord, coordbound):
        
            if coord > 0 and coord < coordbound:
                # if the coordinate is within boundaries:
                return float(coord)
                    # return the original coordinate.
            else:
                return float(0)
                    # reset the coordinate to the beginning of the space frame.
        
        def Data_storage(self, dict1, tcall):
            self.timestore = self.timestep - tcall
                # Time (in time steps) for which the call has been stored into 
                # the dictionary.
            
            if self.timestore > self.max_timestore/self.d_t:
                # if the storing time is longer than it should: 
                dict1[int(self.ID)].pop(tcall, None)
                    # erase it from the memory / dictionary
        
            return dict1

        def Propagation(self, dict1, tmstp):
            self.propdist = float(self.speedsound * self.timestore * env.d_t)
                # propagation distance at timestep according to the time when 
                # the call was emitted and the speed of sound.
            dict1[int(self.ID)][int(tmstp)]['propdist'] = self.propdist
                # update the propagation distance in dict1
        
            return dict1
    
        def In_ring(self, center_x, center_y, radius, x, y):
            dist = (x - center_x) ** 2 + (y - center_y) ** 2
                # distance between an agent and a call's source
            return dist > radius and dist < radius - self.ring_width
                # boolean. Is dist within the distance travelled by the call 
                # between the beginning (radius) and the end of the call 
                # (radius - width)?
        
        def Hearing_test(self, ag_id, tcall, hear_array):
            self.ringtest = self.In_ring(env.callsources[int(ag_id)][int(tcall)]['xsource'],
                                           env.callsources[int(ag_id)][int(tcall)]['ysource'],
                                           env.callsources[int(ag_id)][int(tcall)]['propdist'],
                                           all_bats[int(self.ID)].xhistory[int(self.timestep)],
                                           all_bats[int(self.ID)].yhistory[int(self.timestep)])
                # boolean. Is dist within the distance travelled by the call 
                # between the beginning (radius) and the end of the call 
                # (radius - width)? 
                                           
            if self.ringtest:
                # if focal agent 'self.ID' can hear any call from agent 
                # 'identity' (NB identity = ID is possible, in which case the 
                # bat hears itself):
                hear_array = np.append(hear_array,[self.timestep, ag_id, tcall])
                    # store the ID of the calling bat, the time at which the 
                    # bat called & the time when the call has been heard.
        
            return hear_array

def Dict_update(dict1, dict2):
    # Updates a dictionary without over-writing the keys already stored.
    
    for key in dict2:
        
        if key in dict1:
            dict1[int(key)].update(dict2[int(key)])
        else:
            dict1[int(key)] = dict2[int(key)]
            
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

POPSIZE = 2
    # number of agents to run in one simulation
BOXSIZE = [300,300]
    # spatial boundaries (in meters) within which the agents can move
    # here, the bats can move within an area of 9 ha
DELTA_T = 2
    # time resolution, i.e. time (in ms) represented in 1 simulation time step.
    # Real duration = TIMEFACTOR * simulation duration.
    # allows to keep a sensible ratio & time resolution between pseudo real time and 
    # number of iterations
SIMDURATION = 900

MOVDIRECTION = 0
    # flight direction of the agents
FLIGHTSPEED = 5.5   
    # bats' flight speed in m/s. 5.5 m/s corresponds to a slow bat
    # Hayward & Davis (1964), Winter (1999).
INTER_PULSE_INTERVAL = 50 
    # IPI (ms).  
MAX_HEAR_DIST = 100 
    # maximum distance (in meters) at which a call can be heared
    # ideally, should implement hearing threshold instead e.g. 0 or 20 dB peSPL

# Set the simulation environment with Launcher, according to the parameters given 
# above, and stores it into an object called env
# Run the functions contained in the Launcher class
env = Launcher(POPSIZE, BOXSIZE, DELTA_T, SIMDURATION)
env.Identification()
env.Positions()
env.Timeline()

all_bats = {}
    # creates empty dictionaries for every bat instance to be stored
#all_x = {key:{0:0} for key in list(env.all_ID)}
#all_y = {key:{0:0} for key in list(env.all_ID)}

for ID in env.all_ID:
    all_bats[int(ID)] = env.Bat_Jamming_00(ID, MOVDIRECTION,FLIGHTSPEED,INTER_PULSE_INTERVAL, MAX_HEAR_DIST)
        # stores all instances of the class Bat_Jamming_00 within the bat population
    all_bats[int(ID)].x = env.all_initpos[int(ID)][0]
        # initial x coordinate for each instance, taken from env
    all_bats[int(ID)].xhistory = [all_bats[int(ID)].x]
        # stores it in xhistory
    all_bats[int(ID)].y = env.all_initpos[int(ID)][1]
        # initial y coordinate for each instance, taken from env
    all_bats[int(ID)].yhistory = [all_bats[int(ID)].y]
        # stores it in yhistory
    all_bats[int(ID)].firstcall = 0    
        # time step for initiating the first call 
        # set to zero at the moment, for lack of a better idea
    # all_bats[int(ID)].callduration = 3
        # sets the duration of a bat's calls
        # not necessary for now, but might become useful later
    all_bats[int(ID)].callshistory = np.empty([env.simduration,1], dtype = int)
        # creates an empty list to store call times records
    # all_x[int(ID)][0] = all_bats[int(ID)].xhistory[0]
    # all_y[int(ID)][0] = all_bats[int(ID)].yhistory[0]

for timestep in range(env.simduration):
    for ID in env.all_ID:
        all_bats[int(ID)].timestep = int(timestep)
            # indicates the current time step for each instance 
        all_bats[int(ID)].Movement()
            # makes the instance move
        all_bats[int(ID)].Calling()
            # makes the instance call  
        all_bats[int(ID)].Hearing()
            # makes the instance hear
        plt.plot(all_bats[int(ID)].xhistory, all_bats[int(ID)].yhistory, marker = '^')
            # plot all instances movements over time