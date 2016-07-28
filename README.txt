###########################################################
###													 	###
### 		Documentation for sensory_jamming.py		###
###													 	###
###########################################################

###----------LAUNCHER----------###
### Class implementation to initiate a simulatory environment.
### Inputs: popsize, allpos, boxsize, tres, simduration
### popsize: integer. Size of the bat population / Number of agents to run
### allinitpos: numpy array, dimensions popsize * 2. x & y positions of all bats 
### 	in the population at time 0.
### boxsize: list of 2 integers. Area (rectangle) within which the agents can move, 
### 	defined by the lengths of edges in meters.
### tres: float. Time resolution (s), i.e. amount of time (in s) expressed in
### 	1 simulation time step.
### simduration: integer. Total number of iterations over which to run in the simulation. 

### Once initiated, a Launcher object also contains:
### realduration: duration of the simulation in seconds.
### allID: empty array for integers, dimensions 1*popsize. Will contain the 
### 	identification numbers of all the agents to be run in the simulation.
### callsources: empty dictionary. Will contain the sources of every calls, emitted by 
### 	every agent throughout the simulation.

#----------IDENTIFICATION----------#
# function which assigns a unique serial number to every agent in the population
# Inputs: popsize & allID, taken from the method __init__

#----------OUTPUTS----------#
# allID: updated allID list with all the agents' identification numbers.
#----------End of Documentation for Identification----------#

###----------BAT_JAMMING_00----------###
### Class implementation nested into Launcher class implementation.
### By convention, Launcher class implementations are named env, therefore in
### order to call parameters stored in the Launcher, Bat_Jamming_00 class objects
### refer to them as env.[parameter]
### Makes an agent in a simulation move, call and hear "independently" from the others
### Inputs: ID, movdirection, flightspeed, IPI, max_hear_dist
### ID: integer. Identification number of the focal agent
### movdirection: float. Direction (in radians) towards which the agent is flying
### ### Must be between -pi and pi. Angle relative to the x-axis
### flightspeed: float. Speed, in m/s, at which the agent flies
### IPI: float. Inter-pulse interval, i.e. time interval between each call emitted
### ### by the agent.
### max_hear_dist: float. Distance maximal that a sound can travel before its intensity
### ### passes below the hearing threshold of the agent

### Once initiated, a Bat_Jamming_00 object also contains: 
### boxsize: list of 2 elements, taken from env
### d_t: float, taken from env
### stepsize: float. Distance in meters covered by the agent in 1 timestep
### speedsound: float. Speed of sound at sea level i.e. 340.29
### max_timestore: float. Time maximal, in seconds, during which a sound can travel,
### ### before its intensity passes below the hearing threshold of the agent.
### ring_width: float. Width in meters of the acoustic ring / "doughnut", i.e. spatial 
### ### distribution of a propagating call. Set to 1 meters.
### hearhistory: empty array, dimensions 3*unknown. Will store the identification number 
### ### of the agent who's call has been heard by the focal bat; the time at which the 
### call was emitted, and the time at which it has been heard by the focal agent.

#----------MOVEMENT----------#
# Function that makes the focal agent move over 1 time step
# Inputs: stepsize & movdirection taken from the method __init__

#----------OUTPUTS----------#
# newx: float. New abscissa of the focal agent after 1 time step is spent
# newy: float. New ordinate of the focal agent after 1 time step is spent
# x: float. Updated current abscissa, checked for being within the spatial boundaries
# y: float. Updated current ordinate, checked for being within the spatial boundaries
# xhistory: list of n=simduration floats. Will contain the abscissae of the focal agent
# 	for every time step in the simulation.
# yhistory: list of n=simduration floats. Will contain the ordinates of the focal agent
# 	for every time step in the simulation.
#----------End of Documentation for Movement----------#

#----------CALLING----------#
# Function that makes the focal agent call over 1 time step; stores the call and its 
# sources information (identification number of the agent, position when call is emitted);
# updates the propagation information of every other call; 
# and clears the history of too old calls.
# Inputs: timestep, callstarttime, callshistory. IPI, xhistory & 
# 	yhistory are taken from the __init__ method;
# 	callsources is taken from env.
# timestep: integer. Current time step in the simulation
# callstarttime: float. Time of first call of the agent
# callshistory: numpy array of booleans, dimensions 1*simduration. The   

#----------OUTPUTS----------#
# calltest: float. Number of calls emitted by the focal agent since the 
# beginning, according to the current time and the bat's own IPI
# callshistory: updated callshistory list with the new call if one has just been emitted.
# callsources: updated callsources in env, with the new calls emitted by the focal agent
# 	 and without the outdated / too faint calls.
#----------End of Documentation for Calling----------#

#----------Hearing----------#
# Function that makes the focal agent listen to the soundscape & register the sounds it hears.
# Inputs: all_ID, callsources taken from env, hearhistory taken from __init__ method

#----------OUTPUTS----------#
# hearhistory: Updated hearhistory numpy array with information (time, origin) on the new sound heard, if any
#----------End of documentation for Hearing----------#

#----------Boundaries---------#
# Function that checks if the coordinates of an individual are within the defined boundaries
# Inputs: coord, coordbound

#----------OUTPUTS----------#
coord: Updated coordinate or original coordinate, depending whether it was orginally within boundaries
#----------End of Documentation for Boundaries----------#

#----------DATA_STORAGE----------#
# Function that erases call information when the call is too old to be heard anymore
# Inputs: dict1, tcall
# dict1: dictionary. Contains the call to be checked for
# tcall: time at which the call was emitted

#----------OUTPUTS----------#
# dict1: Updated dict1 dictionary, where the call has been erased if too old
#----------End of Documentation for Data_storage----------#

#----------Propagation----------#
# Function that increments the distance of propagation of a call with time
# Inputs:  dict1, tmstp
# dict1: dictionary. Contains the distance information to be updated
# tmstp: time at which the call was emitted

#----------OUTPUTS----------#
# dict1: Updated dict1 dictionary, with the incremented value of propagation distance
#----------End of Documentation for Propagation----------#

#----------IN_RING----------#
# Function that tests whether a certain point (agent) is on the ring of sound.
# Inputs: center_x, center_y, radius, x, y
# callsource_x: float. Abscissa of the source of the call
# callsource_y: float. Ordinate of the source of the call
# dcfs: float. Stands for Distance of Call from Source. Distance of the propagated call from its source 
# x: float. Abscissa of the object / point to be checked for
# y: float. Ordinate of the object / point to be checked for

#----------OUTPUTS----------#
# dist: distance between an agent and a call source
# This functions returns a boolean. If TRUE, the agent is on the ring
#----------End of Documentation for In_ring----------#

#----------HEARING_TEST----------#
# Function that stores information (time, origin) of a call if heard by the focal agent
# Inputs: ag_id, tcall, hear_array
# ag_id: integer. Identification number of any agent in the simulation7
# tcall: float. Time at which the agent (referred to as ag_id) called
# hear_array: numpy array. Contains the call information (time, origin) of calls heard by the focal agent

#----------OUTPUTS----------#
# ringtest: boolean. TRUE if the focal agent heard the call emitted at tcall by ag_id
# 	FALSE else.
# hear_array: Updated hear_array numpy array with the call information (time, origin) if heard.
#----------End of Documentation for Hearing_test----------#
###----------End of Documentation for Bat_Jamming_00----------###
###----------End of Documentation for Launcher----------###

#----------Dict_update----------#
# Function for updating a dictionary without over-writing the keys already 
# stored in the dictionary.
# Inputs: dict1, dict2
# dict1: original dictionary to be updated.
# dict2: dictionary to add to dict1. NB: has to be of the same format as dict1.

#----------OUTPUTS----------#
# dict1: updated dictionary (new key is added to the already existing keys in the 
#   dictionary, instead of replacing them) 
#----------End of documentation for DICT_UPDATE----------#