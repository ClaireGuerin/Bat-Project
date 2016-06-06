#########################################################################################
### 																				  ###
### Documentation for bat_class program (need to change the name of the programme btw ###
###																					  ###
#########################################################################################

###----------LAUNCHER----------###
### Class implementation to launch the simulation environment.
### By convention, an object of class Launcher() is stored under the name env (stands for environment)
### Initiating inputs: popsize, boxsize, simduration, realduration, all_ID, all_initpos, timecount, timeclock, timeresolution, callsources
### popsize: positive integer. number of agents to simulate i.e. size of the bat group.
### boxsize: list of 2 elements. Cartesian space within which the agents can move
### simduration: positive integer. Duration of simulation 
### ### sim stands for simulation throughout code
### realduration: positive integer. "Real" duration the simulation is supposed to reflect. 
### all_ID: empty list of dimensions 1x(n=popsize elements). 
### all_initpos: empty array of n=popsize lists of 2 elements each. 
### ### init & pos respectively stand for initial and position throughout the code
### timecount: positive integer. Timing of simulation iteration.
### timeclock: array of dimension 1 x simduration. 
### timeresolution: positive float. Resolution of the time simulated.
### callsources: dictionary. Keeps track of the coordinates of all calls emitted by 
###     every agent, at each timestep, as well as the updated propagation distance from the 
###     source at the current time

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

###----------BAT_JAMMING_00----------###
### Class implementation nested into Launcher class for agents jamming:
### - within defined boundaries;
### - in a straight direction;
### - while calling at regular IPI;
### - only timing of sound taken into account;
### - without echoes.
### Initiating Inputs: ID, movdirection, stepsize, IPI, boxsize, simduration, initpos, x, y, callstarttime, xhistory, yhistory, callshistory, hearhistory
### ID: positive integer. ID of the agent.
### movdirection: integer within [-pi;pi]. Angle (in RADIANS) of the agent's direction along the plane. 
### ### mov stands for movement throughout the code.
### stepsize: positive integer. Distance covered by the agent in a single iteration. 
### ### iter & dist respectively stand for iteration & distance throughout the code.
### IPI: inter-pulse interval
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
### hearhistory: numpy array of dimensions n * 4. Reports when an agent heard a sound.
### ### Data provided is, in order: 
### ### - ID of the hearing agent
### ### - time of hearing
### ### - ID of the calling agent = source identification
### ### - time at which the source called

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
### Class implementation for agents hearing:
### - calls from the others and themselves;
### - considering sound propagates in a ring shape, towards all directions
### Initiating inputs: ID, timestep, max_timestore
### ID: integer. Focal bat/agent identification number for which sound hearing is recorded
### timestep: integer. Simulation time step at which the bat is hearing
### max_timestore: integer. Maximum time for storing the data, in simulation steps.
###     Corresponds to the time after which the intensity of the call is considered too low to be heard anymore.
###     Ideally, should implement hearing threshold instead e.g. 20 (pe)SPL

#----------TOA----------#
# Function recording the calls heard by a specific agent, at a soecific time.
# Inputs: ID, timestep, max_timestore, callsources, hearhistory taken from __init__ method

#----------OUTPUTS----------#
# callsources: dictionary. Updated record of all calls sources
# callshearing: numpy array. Updated record of all heard calls
#----------End of documentation for TOA----------#

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
###----------end of documentation for BAT_JAMMING_00----------###