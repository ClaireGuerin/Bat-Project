# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 14:31:31 2016

@author: Claire
"""

from __future__ import division
# import random as rd
import math  as m
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pth

# import pylab

class Launcher:
    # Class implementation to initiate simulatory environment.

    def __init__(self, popsize, boxsize, t_res, simduration):
        self.popsize = int(popsize)
            # integer. Size of the bat population / Number of agents 
        self.t_res = float(t_res)
            # float. Time resolution, i.e. how much time (in s) 
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
        self.realduration = float(self.simduration * self.t_res)
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
        # function which assigns a unique serial number to each agent in the population        
        for agent in range(self.popsize):
            # for each agent in the bat population:
            self.all_ID[int(agent)] = agent
                # fill-in a list of all agents ID 
                # from 0 to n = (population size - 1)
            
    def Positions(self):
        
        for ID in self.all_ID:
            # for each agent in the bat population:
           self.agent_initx = raw_input("Please enter bat #%s x-position in [0-%s]: " % (ID, self.boxsize[0]))
               # ask user for a each agent's coordinate on the x-axis, 
               # within the environment boundaries
           self.agent_inity = raw_input("Please enter bat #%s y-position in [0-%s]: " % (ID, self.boxsize[1]))
               # ask user for a each agent's coordinate on the y-axis, 
               # within the environment boundaries
           self.all_initpos[int(ID),0] = self.agent_initx
               # store the x-coordinate in an array
           self.all_initpos[int(ID),1] = self.agent_inity
               # store the y-coordinate in an array
    
    def Timeline(self):
            
        for iteration in range(self.simduration):
            # for each time step in the simulation:
            self.realtime += float(self.t_res)
                # increment timecount with time resolution. 
            self.timeclock[int(iteration)] = self.realtime
                # store the new, incremented time count in timeclock.

    class Bat_Jamming_00:
        # Class implementation nested into Launcher class implementation.
        # Sets an agent's movement, calling and hearing rules.
    
        def __init__(self, ID, movdirection, flightspeed, IPI, max_hear_dist, callduration):
            self.ID = int(ID)
                # integer. Identification number of the focal agent.
            self.boxsize = env.boxsize
                # takes boxsize from env (class Launcher).
            self.t_res = env.t_res
                # takes t_res from env (class Launcher).
            self.movdirection = float(movdirection)
                # float. Direction (in radians) towards which the agent flies.
                # Must be between -pi and pi (asserted later).
            self.stepsize = float(flightspeed * self.t_res) 
                # float. Distance in meters covered by the agent in 1 time step
            self.IPI = np.around(float(IPI / env.t_res),0)
                # int. Agent's inter-pulse interval, converted in time steps.
            self.speedsound = 340.29
                # float. Speed of sound at sea level in m/s.
            self.max_timestore = float(max_hear_dist / self.speedsound)
                # float. Time maximal, in seconds, during which a sound can 
                # travel, before its intensity passes below the hearing 
                # threshold of the agent.
            self.callduration = callduration
            self.ring_width = float(self.callduration * self.speedsound) 
                # float. difference in radii of the two concentric circles 
                # (in 2D) which form the start and the end of the bat call.
            self.hearhistory_t = []
            self.hearhistory_i = []
            self.hearhistory_c = []
                # empty arrays, length: Number of sounds heard. Will store: 
                # - _i: ID of the source agent's call
                # - _t: time at which the call was emitted 
                # - _c: time at which it was heard by the focal agent        
            assert self.movdirection <= m.pi and self.movdirection >= -(m.pi), "'movdirection' must be in radians & comprised between -pi & pi."
                # returns an error message if movdirection is not within [-pi;pi].
            assert self.IPI > self.t_res, "Inter-pulse interval has to be larger than the time resolution"            
            
        def Movement(self):
            self.newx = float(self.x + self.stepsize * m.cos(self.movdirection))
                # calculate the new x coordinate according to:
                # - the distance travelled over 1 time step.
                # - the direction of the movement
            self.newy = float(self.y + self.stepsize * m.sin(self.movdirection))
                # calculate the new y coordinate according to:
                # - the distance travelled over 1 time step.
                # - the direction of the movement
        
            self.x = self.newx
                # no closed boundaries conditions on x axis.
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
                    
                    self.all_calltimes = env.callsources[int(identity)].keys()
                    self.all_calltimes.sort()
                    
                    for calltime in self.all_calltimes:
                        # for each time step at which the agent previously called:
                        self.Hearing_test(identity, calltime)
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
            
            if self.timestore > self.max_timestore/self.t_res:
                # if the storing time is longer than it should: 
                dict1[int(self.ID)].pop(tcall, None)
                    # erase it from the memory / dictionary
        
            return dict1

        def Propagation(self, dict1, tmstp):
            self.propdist = float(self.speedsound * (self.timestore + 1) * env.t_res)
                # propagation distance at timestep according to the time when 
                # the call was emitted and the speed of sound.
            dict1[int(self.ID)][int(tmstp)]['propdist'] = self.propdist
                # update the propagation distance in dict1
        
            return dict1
    
        def In_ring(self, callsource_x, callsource_y, dcfs, x, y):
            # function which tests if a bat is within the 'ring of sound' of a call
            dist = float(m.sqrt((x - callsource_x) ** 2 + (y - callsource_y) ** 2))
                # distance between an agent and a call's source
            return dist <= dcfs and dist >= dcfs - self.ring_width
                # boolean. Is dist within the distance travelled by the call 
                # between the beginning (dcfs) and the end of the call 
                # (dcfs - width)?
        
        def Hearing_test(self, ag_id, tcall):
            callcentre_x = env.callsources[int(ag_id)][int(tcall)]['xsource']
            callcentre_y = env.callsources[int(ag_id)][int(tcall)]['ysource']
            beam_radius = env.callsources[int(ag_id)][int(tcall)]['propdist']
            agent_xpos = all_bats[int(self.ID)].xhistory[int(self.timestep)]
            agent_ypos = all_bats[int(self.ID)].yhistory[int(self.timestep)]
            
            
            self.ringtest = self.In_ring(callcentre_x, callcentre_y, beam_radius, agent_xpos, agent_ypos)
                # boolean. Is dist within the distance travelled by the call 
                # between the beginning (radius) and the end of the call 
                # (radius - width)? 
                                           
            if self.ringtest:
                # if focal agent 'self.ID' can hear any call from agent 
                # 'identity' (NB identity = ID is possible, in which case the 
                # bat hears itself):
                self.hearhistory_t = np.append(self.hearhistory_t, self.timestep)
                    # store the time when the call has been heard                
                self.hearhistory_c = np.append(self.hearhistory_c, tcall)
                    # store the time at which the call has been emitted
                self.hearhistory_i = np.append(self.hearhistory_i, ag_id)
                    # store the ID of the calling bat 
        
            return self.hearhistory_t, self.hearhistory_c, self.hearhistory_i

def Dict_update(dict1, dict2):
    # Updates a dictionary without over-writing the keys already stored.
    
    for key in dict2:
        
        if key in dict1:
            dict1[int(key)].update(dict2[int(key)])
        else:
            dict1[int(key)] = dict2[int(key)]
            
    return dict1

### SIMULATIONS ###
        
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
BOXSIZE = [200,200]
    # spatial boundaries (in meters) within which the agents can move
    # here, the bats can move within an area of 9 ha
TIME_RESOLUTION = 0.0360
    # time resolution, i.e. time (in s) represented in 1 simulation time step.
    # Real duration = TIMEFACTOR * simulation duration.
    # allows to keep a sensible ratio & time resolution between pseudo real time and 
    # number of iterations
SIMDURATION = 20

MOVDIRECTION = 0
    # flight direction of the agents
FLIGHTSPEED = 5.5   
    # bats' flight speed in m/s. 5.5 m/s corresponds to a slow bat
    # Hayward & Davis (1964), Winter (1999).
INTER_PULSE_INTERVAL = 0.05 
    # IPI (s).  
MAX_HEAR_DIST = 300 
    # maximum distance (in meters) at which a call can be heared
    # ideally, should implement hearing threshold instead e.g. 0 or 20 dB peSPL
CALL_DURATION = 0.01
    # duration of the call (in seconds)

# Set the simulation environment with Launcher, according to the parameters given 
# above, and stores it into an object called env
# Run the functions contained in the Launcher class
env = Launcher(POPSIZE, BOXSIZE, TIME_RESOLUTION, SIMDURATION)
env.Identification()
env.Positions()
env.Timeline()

all_bats = {}
    # creates empty dictionaries for every bat instance to be stored
#all_x = {key:{0:0} for key in list(env.all_ID)}
#all_y = {key:{0:0} for key in list(env.all_ID)}

for ID in env.all_ID:
    all_bats[int(ID)] = env.Bat_Jamming_00(ID, MOVDIRECTION,FLIGHTSPEED,INTER_PULSE_INTERVAL, MAX_HEAR_DIST, CALL_DURATION)
        # stores all instances of the class Bat_Jamming_00 within the bat population
    all_bats[int(ID)].x = env.all_initpos[int(ID)][0]
        # initial x coordinate for each instance, taken from env
    all_bats[int(ID)].xhistory = [all_bats[int(ID)].x]
        # stores it in xhistory
    all_bats[int(ID)].y = env.all_initpos[int(ID)][1]
        # initial y coordinate for each instance, taken from env
    all_bats[int(ID)].yhistory = [all_bats[int(ID)].y]
        # stores it in yhistory
    all_bats[int(ID)].firstcall = int(0)   
        # time step for initiating the first call 
        # set to zero at the moment, for lack of a better idea
    # all_bats[int(ID)].callduration = 3
        # sets the duration of a bat's calls
        # not necessary for now, but might become useful later
    all_bats[int(ID)].callshistory = np.empty([env.simduration,1], dtype = int)
        # creates an empty list to store call times records 
    
for timestep in range(env.simduration):
    
    for ID in env.all_ID:
        all_bats[int(ID)].timestep = int(timestep)
            # indicates the current time step for each instance 
        all_bats[int(ID)].Calling()
            # makes the instance call
        
        for n in env.callsources[int(ID)].keys():
            callcentre_x = env.callsources[int(ID)][n]['xsource']
            callcentre_y = env.callsources[int(ID)][n]['ysource']
            beam_radius = env.callsources[int(ID)][n]['propdist']
            pth.Circle([callcentre_x,callcentre_y], beam_radius, fill = False)
    
    for ID in env.all_ID:
        all_bats[int(ID)].timestep = int(timestep)
            # indicates the current time step for each instance
        all_bats[int(ID)].Hearing()
            # makes the instance hear
    
    for ID in env.all_ID:
        all_bats[int(ID)].timestep = int(timestep)
            # indicates the current time step for each instance
        all_bats[int(ID)].Movement()
            # makes the instance move
        plt.plot(all_bats[int(ID)].xhistory, all_bats[int(ID)].yhistory, marker = '^')
            # plot all instances movements over time

### EXPORTING DATA ###

filenamesH = []
filenamesM = []

for ID in env.all_ID:
    tmstp = all_bats[ID].hearhistory_t
    tcall = all_bats[ID].hearhistory_c
    idbat = all_bats[ID].hearhistory_i
    all_bats[ID].hearhistory = {"t": tmstp, "c": tcall, "i": idbat}
    
    xtrack = all_bats[ID].xhistory
    ytrack = all_bats[ID].yhistory
    all_bats[ID].movhistory = {"x": xtrack, "y": ytrack}
    
    for data_type in all_bats[ID].hearhistory.keys():
        filenamesH.append("C:\Users\dlewanzik\Documents\Bat-Project\Res\Hearing\%s_hearhistory_%s.txt" % (str(ID), data_type))
    
    for coordinate in all_bats[ID].movhistory.keys():
        filenamesM.append("C:\Users\dlewanzik\Documents\Bat-Project\Res\Moving\%s_movhistory_%s.txt" % (str(ID), coordinate))
        
    with open("C:\Users\dlewanzik\Documents\Bat-Project\Res\Calling\%s_callshistory.txt" % ID, "w") as fp3:
        for value in all_bats[ID].callshistory:
            fp3.writelines('%s\n' % value[0])
    fp3.close()

for fname in filenamesH:
    with open("%s" % fname, "w") as fp1:
        end_id = [n for n in xrange(len(fname)) if fname.find('_', n) == n][0]
        for value in all_bats[int(fname[53:end_id])].hearhistory[fname[-5]]:
            fp1.writelines('%s\n' % value)
        fp1.close()

        
for fname in filenamesM:
    with open("%s" % fname, "w") as fp2:
        end_id = [n for n in xrange(len(fname)) if fname.find('_', n) == n][0]
        for value in all_bats[int(fname[52:end_id])].movhistory[fname[-5]]:
            fp2.writelines('%s\n' % value)
    fp2.close()

### TESTS ###

x1=all_bats[0].xhistory[0]
x2=all_bats[1].xhistory[0]

y1=all_bats[0].yhistory[0]
y2=all_bats[1].yhistory[0]

b_dist = m.sqrt((x1-x2)**2+(y1-y2)**2)

np.where(all_bats[0].callshistory == 1)
np.where(all_bats[1].callshistory == 1)

all_bats[0].hearhistory_t
all_bats[0].hearhistory_i
all_bats[0].hearhistory_c

all_bats[1].hearhistory_t
all_bats[1].hearhistory_i
all_bats[1].hearhistory_c

x11=all_bats[0].xhistory[0:5]
x21=all_bats[0].xhistory[1:6]
x21-x11

x12=all_bats[1].xhistory[0:5]
x22=all_bats[1].xhistory[1:6]
x22-x12


radius = 1
center_x = 1
center_y = 1

pth.Circle([center_x,center_y], radius)

