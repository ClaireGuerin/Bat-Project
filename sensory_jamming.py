# -*- coding: utf-8 -*-
'''
Created on Mon Apr 04 14:31:31 2016
@author: Claire
'''

# %reset
# Reset the console and clear it from all previous variables that might have 
# been stored. Type 'y' to confirm resetting.

from __future__ import division
import random as rd
import math  as m
import numpy as np
import matplotlib.pyplot as plt
import itertools as it
from scipy.optimize import minimize
import itertools
import os, errno


class Launcher:
# Class implementation to initiate simulatory environment.
    
    def __init__(self, tres, simduration):
         
        self.tres = float(tres)
        # float. Time resolution, i.e. how much time (in s) is expressed in
        #  1 simulation time step.
        self.simduration = int(simduration)
        # integer. Simulation duration, i.e. total number of time steps to be 
        # run in the simulation. 
        self.realduration = float(self.simduration * self.tres)
        # float. Duration of the simulation (in s).
        self.callsources = {}
        # empty dictionary. Will contain the sources of each calls, 
        # emitted by every agent throughout the simulation.
        self.speedsound = 340.29
        # float. Speed of sound at sea level in m/s.

    def Square_lattice(self, lowvertex, axIID, Nedge):
        
        self.popsize = int(Nedge ** 2)
        # integer. Size of the bat population / Number of agents.
        self.allID = range(self.popsize)
        # list of floats, dimensions 1*popsize. IDs a unique ID corresponding to the agent.
        x = np.linspace(lowvertex[0], lowvertex[0]+(Nedge-1)*axIID, num=Nedge)
        y = np.linspace(lowvertex[1], lowvertex[1]+(Nedge-1)*axIID, num=Nedge)
        
        #  the initial positions are not exact numbers -- but more like float approxmns - next thing to be figured out !               
        
        bothedges = [x,y]
        self.allinitpos = list(it.product(*bothedges))
        # list of tuples, dimensions popsize*2. 
        # x & y positions of all bats in the population at time 0.

    class Bat_Jamming_00:
    # Class implementation nested into Launcher class implementation.
    # Sets the individual rules for the agents to move, call & hear  
    # during the simulation.
    
        def __init__(self, ID, movangle, flightspeed, IPI, maxheardist, callduration):
            
            self.ID = int(ID)
            # integer. Serial (identification) number of the focal agent.
            self.tres = env.tres
            # take tres from env (class Launcher).
            self.movangle = m.radians(float(movangle))
            # float. Movement angle, i.e. direction (in rad) towards which the 
            # agent flies. Value between 2*pi and pi, as asserted below.            
            self.stepsize = float(flightspeed * self.tres) 
            # float. Distance in meters covered by the agent in 1 time step
            #self.callduration = callduration
            self.callduration= np.around(float(callduration/ env.tres), 0)
            # float. Duration (in s) of each call. 
            #self.IPI = np.around(float((self.callduration/dutycycle - self.callduration) / env.tres),0)
            self.IPI = np.around(float(IPI / env.tres), 0)
            # integer. Inter-pulse interval of the agent, converted in time steps.
            self.maxtimestore = float(maxheardist / env.speedsound)
            # float. Maximum time for storing a sound, deduced from the time 
            # maximal, in seconds, during which a sound can travel, before its 
            # intensity passes below the hearing threshold of the agent.
            
            self.ringwidth = float(self.callduration*self.tres * env.speedsound) 
            # float. Difference in radii of the two concentric circles 
            # (in 2D) which form the start and the end of the bat call.
            self.firstcall = rd.randint(0,self.IPI+self.callduration)
            # time step for initiating the first call
            self.hearhistory_t = []
            self.hearhistory_i = []
            self.hearhistory_c = []
            # empty arrays, length: Number of sounds heard. Will store: 
            # - _i: ID of the agent's call who emitted the call that is heard;
            # - _c: time at which the call was emitted;
            # - _t: time at which it was heard by the focal agent.

            assert self.movangle >= m.radians(0) and self.movangle < m.radians(360), "Enter movement angle between 0 and 360 degrees"
            #assert self.IPI > self.tres, 'Inter-pulse interval must be larger than the time resolution'            
            
        def Movement(self):
            
            self.newx = float(self.x + self.stepsize * m.cos(self.movangle))
            # calculate the new x coordinate according to:
            # - the distance travelled over 1 time step.
            # - the direction of the movement
            self.newy = float(self.y + self.stepsize * m.sin(self.movangle))
            # calculate the new y coordinate according to:
            # - the distance travelled over 1 time step.
            # - the direction of the movement
        
            self.x = self.newx
            # update the current x coordinate of the agent
            self.y = self.newy
            # update the current y coordinate of the agent
        
            self.xhistory = np.append(self.xhistory, self.x)
            # store the updated x coordinate in xhistory
            self.yhistory = np.append(self.yhistory, self.y)
            # store the updated y coordinate in yhistory
            
        def Calling(self):
            
            self.calltest = float(self.timestep-self.firstcall)/float(self.IPI + self.callduration)
            # theoretical number of calls since the 1st call:
            # time since 1st call = current time - time of 1st call/IPI
            
            if self.calltest%1 == 0:
            # if the theoretical number of calls since the 1st call is a 
            # divisible by the call cycle duration then initiate a call
            
                self.callshistory[int(self.timestep)] = 1
                # account for the initiation of a new call at this time step 
                # by adding a 1 in the calls' history.
                Dict_update(env.callsources, {self.ID:{self.timestep:{'xsource': self.xhistory[self.timestep], 'ysource': self.yhistory[self.timestep], 'propdist':0.0}}})
                # store the position of the bat at the time of calling into a 
                # dictionary of the form: {bat:{time of calling:[x,y,propdist]}}
                
            else:
                self.callshistory[int(self.timestep)] = 0
                # account for the absence of call initiation at this time 
                # step by adding a zero in the calls' history.
                
            if self.ID in env.callsources.keys():
            # if the agent has call already stored in the call dictionary:
                self.Sound_update(env.callsources, self.ID)
                # update the information on this call.

        def Hearing(self):
            
            for identity in env.allID:
            # for each agent in the simulation:
                
                if identity in env.callsources.keys():
                # if the agent has a call stored in the call dictionary: 
                    
                    
                    # store the time of emission of the call 
                    self.allcalltimes = env.callsources[identity].keys()
                    self.allcalltimes.sort()
                    
                    
                    for calltime in self.allcalltimes:
                    # for each time step at which any agent has previously called:
                        self.Hearing_test(env.callsources, identity, calltime)
                        # identify and record calls that can be heard by focal agent
        
        # calls_dict is the dictionary which holds all the calls to be propagated in the simulation environment
        
        def Sound_update(self, calls_dict, identity):
            
            allsoundtimes = calls_dict[identity].keys()
            allsoundtimes.sort()
            # emission time of every sound that is recorded in calls_dict for this 
            # individual
            
            for soundtime in allsoundtimes:
            # for each time step at which a sound was emitted:
                self.Data_storage(calls_dict, int(soundtime))
                # erase sounds that are too old to be heared anymore
                
                if soundtime in calls_dict[identity].keys():
                # if the sound is still in the dictionary:
                    self.Propagation(calls_dict, int(soundtime), identity)
                    # update its propagation distance accordingly            
        
        def Data_storage(self, calls_dict, tcall):
            
            self.timestore = self.timestep - tcall
            # Time (in time steps) for which the call has been stored into 
            # the dictionary.
            
            if self.timestore > self.maxtimestore/self.tres:
            # if the storing time is longer than it should: 
                calls_dict[int(self.ID)].pop(tcall, None)
                # erase it from the memory / dictionary
        
            return calls_dict

        def Propagation(self, calls_dict, tmstp, identity):
            
            dictionary = calls_dict[identity][int(tmstp)]
            
            self.backradius = dictionary['propdist']
            # float. Current value for the radius of the sound from its source,
            # that was calculated at the previous timestep.
            self.propdist = float(env.speedsound * (self.timestore + 1) * env.tres)
            # Current propagation distance at timestep according to the time when 
            # the call was emitted, and the speed of sound.
            dictionary['propdist'] = self.propdist
            # update the propagation distance in calls_dict
            
            return calls_dict
    
        def In_ring(self, xcallsource, ycallsource, soundback, soundfront, x, y):
            
            dist = float(m.sqrt((x - xcallsource) ** 2 + (y - ycallsource) ** 2))
            # distance between the agent and the source of the call
            return dist <= soundfront and dist > soundfront - self.ringwidth
            # boolean. Is dist within the distance travelled by the call 
            # between the beginning of the call at t = tres * (timestep - 1) 
            # (soundfront), & the end of the call at t = tres * timestep 
            # (soundback - self.ringwidth)?            
        
        def Hearing_test(self, calls_dict, agID, temission):
            
             # float. x-coordinate of the source of the call emitted by agID at 
            # timestep = temission .
            xcallcentre = calls_dict[int(agID)][int(temission)]['xsource']
            
            # float. y-coordinate of the source of the call emitted by agID at 
            # timestep = temission.    
            ycallcentre = calls_dict[int(agID)][int(temission)]['ysource']
            
            # float. Current propagation distance of the source of the call 
            # emitted by agID at timestep = temission.     
            beamradius = calls_dict[int(agID)][int(temission)]['propdist']
            
            # float. Current value for the radius of the sound from its source,
            # that was calculated at the previous timestep.            
            backradius = beamradius - env.speedsound * env.tres
            
            # float. Current x-coordinate of the focal agent. 
            xposagent = allbats[int(self.ID)].xhistory[int(self.timestep)]
            
            # float. Current x-coordinate of the focal agent.
            yposagent = allbats[int(self.ID)].yhistory[int(self.timestep)]
            
            
            self.ringtest = self.In_ring(xcallcentre, ycallcentre, backradius, beamradius, xposagent, yposagent)
            # boolean. Does the focal agent hear the call emitted by agID
            # at timestep = temission, knowing its current position?
                                           
            if self.ringtest:
                # if the focal agent (self.ID) can hear the call emitted by  
                # agID (NB: agID = ID is possible, in which case the 
                # bat hears itself):
                self.hearhistory_t = np.append(self.hearhistory_t, self.timestep)
                # store the time at which the call has been heard by the focal agent                
                self.hearhistory_c = np.append(self.hearhistory_c, temission)
                # store the time at which the call was emitted
                self.hearhistory_i = np.append(self.hearhistory_i, agID)
                # store the ID of the bat who emitted the call
                
            return self.hearhistory_t, self.hearhistory_c, self.hearhistory_i


def Dict_update(calls_dict, dict2):
# Updates a dictionary without over-writing the keys already stored in it.
    
    for key in dict2:
    
        if key in calls_dict:
            calls_dict[key].update(dict2[key])
        else:
            calls_dict[key] = dict2[key]
            
    return calls_dict

def Min_hear(param,alpha,sourcelevel, hth):
# function to calculate the maximum distance at which a bat call can be heard
    SL = sourcelevel - m.fabs(20 * m.log10(0.1/1)) # source level corrected for intensity at 1m
    hearsqr = (SL - 20 * m.log10(param) + alpha * param - hth) ** 2
    return hearsqr

### SIMULATIONS ###      
### Run a multi-agents simulation and plot agents' movements.

###----------PARAMETERS----------###
### to be determined / chosen beforehand

### GENERAL ### 

### POPULATION_SIZE: integer. Size of the population of bats you want to simulate
### INITIAL_POSITIONS: numpy array, dimensions 2 * population size. Contains the
###     initial x and y coordinates of each agent in the population, in meters (m).
### BOX_SIZE: list of 2 float elements. Space within which the bats are moving,
###     in meters (m).
### TIME_RESOLUTION: float. Time resolution, i.e time (in s) represented over 1 
###     simulation time step.allows to keep a sensible ratio between time in
###     in seconds & time in time steps in the simulation.
###     Time (in s) = TIME_RESOLUTION * time (in simulation time steps).
### SIMULATION_DURATION: duration of the simulation, i.e. number of iterations
###     over which the simulation must be ran.

### INDIVIDUAL ###

### MOVEMENT_ANGLE: float. Angle (in degrees), between the x-axis and the 
###     direction of the movement of the agent. 
###     Only values between 0 & 360°C are accepted.
### FLIGHT_SPEED: float. Flight speed (m/s) of the agent. 
###     5.5 m/s corresponds to a slow bat (Hayward & Davis (1964), Winter (1999)).
### INTER_PULSE_INTERVAL: Inter-pulse interval of the agent, i.e. time interval 
###     between each call initiated.
### DUTY_CYCLE = float. Duty cycle of the agent (ratio IPI:call duration).
### MAXIMUM_HEARING_DISTANCE: float. Maximum distance (m) at which a sound can be
###     heard by the agent (depending on the intensity threshold)

TIME_RESOLUTION = 0.001 
SIMULATION_DURATION =20

CORNER_INDIVIDUAL_POSITION = [1,1]
IID_ON_AXE = 0.9
N_EDGE = 2

MOVEMENT_ANGLE = 0
FLIGHT_SPEED = 0
#DUTY_CYCLE = 0.1228
CALL_DURATION = 0.002
INTER_PULSE_INTERVAL = 0.008
ALPHA = -1.7 # db/m absorption at particular frequency 
SOURCE_LEVEL = 120 # dB SPL, ref 20uPa @10cm 
HEARING_THRESHOLD = -10 #hearing threshold in dB SPL 
MAXIMUM_HEARING_DISTANCE = float(minimize(Min_hear, x0 = 30, args = (ALPHA,SOURCE_LEVEL,HEARING_THRESHOLD))['x'])

rd.seed(96) # initialize the basic random number generator.

# Set the simulation environment with Launcher, according to the given parameters 
# and store it into an object called env.


## -- running one instance of the simulation and getting data outputs :

# PLEASE CREATE A DIRECTORY CALLED 'RES' IN ANY DESIRED LOCATION AND PASTE THE ADDRESS (' C:\\_____\\.......\\Res')

currdir='C:\\Users\\tbeleyur\\Desktop\\Res'

os.chdir(currdir) # change working directory


    
env = Launcher(TIME_RESOLUTION, SIMULATION_DURATION)
env.Square_lattice(CORNER_INDIVIDUAL_POSITION, IID_ON_AXE, N_EDGE)

allbats = {}
# Create an empty dictionary for every bat instance to be stored

# Set the agents with Bat_Jamming_00, according to the given parameters.

for ID in env.allID:

    allbats[int(ID)] = env.Bat_Jamming_00(ID, MOVEMENT_ANGLE,FLIGHT_SPEED,INTER_PULSE_INTERVAL, MAXIMUM_HEARING_DISTANCE, CALL_DURATION)
    # store all instances of the class Bat_Jamming_00 within the bat population
    allbats[int(ID)].x = env.allinitpos[int(ID)][0]
    # initial x coordinate for each instance, taken from env
    allbats[int(ID)].xhistory = [allbats[int(ID)].x]
    # store it in xhistory
    allbats[int(ID)].y = env.allinitpos[int(ID)][1]
    # initial y coordinate for each instance, taken from env
    allbats[int(ID)].yhistory = [allbats[int(ID)].x]
    # store it in yhistory
    allbats[int(ID)].callshistory = np.empty([env.simduration,1], dtype = int)
    # create an empty list to store records of call times

    fig, ax = plt.subplots()
    # create the figure where rings of sound will be drawn & bats positions will be plotted
    colorpanel = plt.get_cmap('nipy_spectral')
    # select a color panel to differentiate individuals

# Run the simulation of the individual based model, & store the results in allbats 

for timestep in range(env.simduration):
    
    for ID in env.allID:
    
        allbats[int(ID)].timestep = int(timestep)
        # current time step for each instance
        allbats[int(ID)].Calling()
        # make the instance call
        #print(timestep)
    
        if ID in env.callsources.keys():
    
            for n in env.callsources[int(ID)].keys():
        
                xcallcentre = env.callsources[int(ID)][n]['xsource']
                ycallcentre = env.callsources[int(ID)][n]['ysource']
                beamradius = env.callsources[int(ID)][n]['propdist']
                ringout = plt.Circle([xcallcentre,ycallcentre], beamradius, color = colorpanel(ID*100), fill = False)
                ax.add_artist(ringout)
                # plot the corresponding ring of sound
        
    for ID in env.allID:
    
        allbats[int(ID)].timestep = int(timestep)
        # current time step for each instance
        allbats[int(ID)].Hearing()
        # make the instance hear

    for ID in env.allID:
    
        allbats[int(ID)].timestep = int(timestep)
        # current time step for each instance
        allbats[int(ID)].Movement()
        # make the instance move
        positions = plt.plot(allbats[int(ID)].xhistory, allbats[int(ID)].yhistory, color = colorpanel(ID*100), marker = '^')
        ax.add_artist(positions[0])
        # plot all instances movements over time

ax.set_xlim(-1,15)
# set the x-limit of the figure
ax.set_ylim(-1,15)
# set the y-limit of the figure
fig.savefig('plot_bats_rings.pdf')
# save the figure
plt.close()
# close the figure



# function which makes a directory if it doesn't already exist
def dirmaker(dirname):
    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname)
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

# creating the necessary directories where data will be saved
dirname = "ipi%s_nedge%s_iidaxe%s" % (str(INTER_PULSE_INTERVAL), str(N_EDGE), str(IID_ON_AXE))
dirmaker(dirname)
dirmaker(os.path.join(dirname,'Moving'))
dirmaker(os.path.join(dirname,'Calling'))
dirmaker(os.path.join(dirname,'Hearing'))




for ID in env.allID:

    allbats[ID].xhistory = np.delete(allbats[ID].xhistory, 0)
    # remove initial x positions, which were use din the algorithm but are 
    # not integrated into movement
    allbats[ID].yhistory = np.delete(allbats[ID].yhistory, 0)
    # remove initial y positions, which were used in the algorithm but are 
    # not integrated into movement





    ### EXPORTING DATA ###

filenamesH = [] # hearing files 
filenamesM = [] #movement files 


for ID in env.allID:

    tmstp = allbats[ID].hearhistory_t
    tcall = allbats[ID].hearhistory_c
    idbat = allbats[ID].hearhistory_i
    allbats[ID].hearhistory = {'t': tmstp, 'c': tcall, 'i': idbat}

    xtrack = allbats[ID].xhistory
    ytrack = allbats[ID].yhistory
    allbats[ID].movhistory = {'x': xtrack, 'y': ytrack}

    for data_type in allbats[ID].hearhistory.keys():
        filenamesH.append('%s\Hearing\%s_hearhistory_%s.txt' % (dirname, str(ID), data_type))

    for coordinate in allbats[ID].movhistory.keys():
        filenamesM.append('%s\Moving\%s_movhistory_%s.txt' % (dirname, str(ID), coordinate))
    
    with open('%s\Calling\%s_callshistory.txt' % (dirname, ID), 'w') as fp3:
        for value in allbats[ID].callshistory:
            fp3.writelines('%s\n' % value[0])
    fp3.close()


# --- function which determines the bat Id in concern when given a filename
import re # import the regular expressions function   

def idfinder(filename):
     # figure out which is the bat ID of the fname in concern
    
    #find all the numbers in the fname string, and take the last number (batID)
    batid=int(re.findall(r'\d+', fname)[-1])
    
    return(batid)
# --end of function 
 
#writing the values of hearhistory for each individual into the file

for fname in filenamesH:
    with open('%s' % fname, 'w') as fp1:
        
        # find the ID number of the bat in concern from the filename        
        batid=idfinder(fname);
        
        #for value in allbats[int(fname[27:end_id])].hearhistory[fname[-5]]:
        for value in allbats[batid].hearhistory[fname[-5]]:
            
            fp1.writelines('%s\n' % value)
            
    fp1.close()


# writing the value of move history into the file for each individual 
    
for fname in filenamesM:
    with open('%s' % fname, 'w') as fp2:
        
        # find the ID number of the bat in concern from the filename
        batid=idfinder(fname);
        
        for value in allbats[batid].movhistory[fname[-5]]:
            
            fp2.writelines('%s\n' % value)
            
    fp2.close()
    

    
