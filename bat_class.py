# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 14:31:31 2016

@author: Claire
"""

import random as rd
import math  as m
import numpy as np
import matplotlib.pyplot as plt
# import pylab

class Launcher:
    
    def __init__(self, popsize, boxsize, timefactor, realduration):
        self.popsize = int(popsize)
        self.timefactor = timefactor
            # factor expressing the difference between simduration and realduration
        self.timecount = 0 
            # set the time at 0.
        self.timeclock = np.empty(self.simduration, dtype = float)
        self.timeclock[0] = float(self.timecount) 
            # record the first time, i.e. 0.
        self.realduration = float(realduration)
            # real duration = time factor * simduration
        self.simduration = float(self.realduration/10)
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
            self.timecount += float(1/self.timefactor)
                # increment timecount with time resolution. 
                # Nothing weird here: Timeline creates a clock in terms of real time
                # As 1 simulation iteration is ran, there is simduration/realduration
                # time spent "for real", i.e. time resolution
                # timecount corresponds to the "actual" time, as opposed to simulated time
            self.timeclock[iteration] = self.timecount
                # store the new, incremented timecount in timeclock.

    class Bat_Jamming_00:
    
        def __init__(self, ID, movdirection, flightspeed, IPI, max_timestore, ring_dev):
            self.ID = int(ID)
            self.boxsize = env.boxsize
                # takes boxsize from env, an object of the class Launcher
            self.simduration = int(env.simduration)
                # takes simduration from env, an object of the class Launcher
            self.movdirection = float(movdirection)
            self.stepsize = float(flightspeed)       
            self.IPI = float(IPI)
            self.max_timestore = float(max_timestore)
            self.ring_dev = int(ring_dev)
            self.callsource = {}
            self.hearhistory = []
        
            assert self.movdirection <= m.pi and self.movdirection >= -(m.pi), "'movdirection' must be in radians & comprised between -pi & pi."
                # returns an error message if movdirection is not within [-pi;pi].
    
        def Movement(self):
            self.newx = self.x + self.flightspeed * float(1/env.timefactor) * m.cos(self.movdirection)
                # calculates the new x coordinate according to:
                # ### the distance travelled over 1 time step.
                # ### the direction of the movement
            self.newy = self.y + self.flightspeed * float(1/env.timefactor) * m.sin(self.movdirection)
                # calculates the new y coordinate according to:
                # ### the distance travelled over 1 time step.
                # ### the direction of the movement
        
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
                # (time since the 1st call = current time - time of the 1st call/IPI)
            
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

        def In_ring(self, center_x, center_y, radius_B, x, y, radius_E):
            dist = (x - center_x) ** 2 + (y - center_y) ** 2
            return dist > radius_B and dist < radius_E
                # the distance between an agent and a call's source has to within 
                # the distance travelled by the call between the beginning (radius_B)
                # and the end of the call (radius_E).
        
        def Propagation(self, dict1, tmstp):
            SPEED_SOUND = 340.29 # speed of sound at sea level in m/s
            self.propdist = float(SPEED_SOUND * self.timestore * env.timefactor)
                # calculate, for a certain call, its propagation distance at timestep 
                # according to the time when the call was emitted and the speed of sound
            dict1[self.ID][tmstp]['propdist'] = self.propdist
                # update the propagation distance in dict1
        
            return dict1
    
        def Hearing_test(self, ag_id, tcall, hear_array, tmstp):
            self.ringtest = self.Min_circle(self.callsource[ag_id][tcall]['xsource'],
                                           self.callsource[ag_id][tcall]['ysource'],
                                           self.callsource[ag_id][tcall]['propdist'],
                                           all_bats[int(ID)].xhistory[tmstp],
                                           all_bats[int(ID)].yhistory[tmstp],
                                           self.callsource[ag_id][tcall]['propdist'] + self.callduration)
                                         
            if self.ringtest:
                # test if the agent 'ID' can hear any call from agent 'identity'
                # including own calls
                hear_array = np.append(hear_array,[self.ID, tmstp, ag_id, tcall])
                    # if so, store the IDs of the hearing bat and the calling bat,
                    # as well as the time at which the bat called and has been heard
        
            return hear_array

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
    # number of agents to run in one simulation
BOXSIZE = [2000,2000]
    # spatial boundaries (in meters) within which the agents can move
TIMEFACTOR = 10
    # Real duration = TIMEFACTOR * simulation duration.
    # the reason why I set it like this is that it allows to keep a sensible ratio
    # & "time resolution" between pseudo real time and number of iterations
REALDURATION = 10000
    # total time in milliseconds for the duration of simulation (# iterations)
MOVDIRECTION = 0
    # flight direction of the agents
FLIGHTSPEED = 5.5   
    # bats' flight speed in m/s. 5.5 m/s corresponds to a slow bat
    # Hayward & Davis (1964), Winter (1999).
INTER_PULSE_INTERVAL = 50 
    # IPI in seconds (s).  
MAX_TIMESTORE = 100 
    # ideally, should implement hearing threshold instead e.g. 0 or 20 dB peSPL
RING_DEV = 50
    # to be changed.

# Set the simulation environment with Launcher, according to the parameters given 
# above, and stores it into an object called env
# Run the functions contained in the Launcher class
env = Launcher(POPSIZE, BOXSIZE, TIMEFACTOR, REALDURATION)
env.Identification()
env.Positions()
env.Timeline()

all_bats = {}
    # creates empty dictionaries for every bat instance to be stored
#all_x = {key:{0:0} for key in list(env.all_ID)}
#all_y = {key:{0:0} for key in list(env.all_ID)}
all_sources = {}

for ID in env.all_ID:
    all_bats[int(ID)] = env.Bat_Jamming_00(int(ID), MOVDIRECTION,FLIGHTSPEED,INTER_PULSE_INTERVAL, MAX_TIMESTORE, RING_DEV)
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
    all_bats[int(ID)].callduration = 3
        # sets the duration of a bat's calls
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