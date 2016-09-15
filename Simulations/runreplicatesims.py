# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 12:57:40 2016

@author: tbeleyur
"""
import sys
import os
import numpy as np
import random as rd
# script which runs multiple replicates of many parameter combinations:
# inputs : 
# results_dir: folder in which ALL the results will be stored
# Nrep : number of replicate simulations to be run per parameter set
# Param_set: list with sublists. [ [parcombi1],[parcombi2],parcombi3..] 

#outputs:
# results in the specificed folder - each parameter combination gets its own labelled folder within each 
# replicate numbered folder (Res#repnum)

# uncomment to see how long the whole script takes (also uncomment last line)
import time
start = time.clock()


# setting random seed :

rd.seed(01) # initialize the basic random number generator - which should make sure the results are compatible across runs and multiple systems! 


#the target folder that you want all the results to be stored in : (folder doesn't need to exist ! ) 
Results_dir='C:\\Users\\qliang\\Documents\\Bat-Project'


# import the script ('oneinst.py') which runs one single instance of the simulation : 
modulelocn='C:\\Users\\qliang\\Documents\\GitHub\\Bat-Project'

sys.path.append ( modulelocn ) # add the location of the module to the search path of python 
from oneinst import onerun  # import the one function from the module 


#the index numbers of the parameters being used in the pcomb lists : 
    
    #0 N_SIDE: integer. the number of bats per side in a square lattice - the total population size will be N_SIDE^2
    #1 TIME_RESOLUTION: float. Time resolution, i.e time (in s) represented over 1 iteration
    #     simulation time step.allows to keep a sensible ratio between time in
    #     in seconds & time in time steps in the simulation.
    #     Time (in s) = TIME_RESOLUTION * time (in simulation time steps).
    #2 SIMULATION_DURATION: integer. duration of the simulation, i.e. number of iterations
    #3     over which the simulation must be ran.
    #4 CORNER_INDIVIDUAL_POSITION : 1x2 list with float numbers, the 'corner' position from which the rest of the lattice is built
    #5 IID_ON_AXE : float, inter individual distance on the axis 
       
    #6 MOVEMENT_ANGLE: float. Angle (in degrees), between the x-axis and the 
    #     direction of the movement of the agent. 
    #     Only values between 0 & 360°C are accepted.
    #7 FLIGHT_SPEED: float. Flight speed (m/s) of the agent.
    #8 CALL_DURATION: float. length of call, in seconds
    #9 INTER_PULSE_INTERVAL: Inter-pulse interval of the agent, i.e. time interval 
    #     between each call initiated.
    #10 HEARING_THRESHOLD: lowest sound pressure level a bat can hear (dB SPL)
    #11 SOURCE_LEVEL : the sound pressure level of a bat call (dB SPL @ 10cm)
    #12 ALPHA : atmospheric absorption of sound at the call frequencies
    #13 speed of sound : float. velocity of sound propagation in m/s
 

 
# display the simulation scenario being run here - please change depending on the variable being varied:
 
SIMULATION_GOAL='to study changes in ipi unconstrained'

print('\n %s \n'%SIMULATION_GOAL)
   


nedge=2
tres=0.001
sim_durn=0 # dummy value which is later changed 
cornerpos=(1,1)
iid_on_axis=0.6
mov_angle=0
flight_speed=4
call_durn=0.002

ipi_combs=0.004

hear_thresh=-10
source_level=120
alpha=-1.41
speedsound=341.5





# create the parameter combinations wanted : - right now it is purposefully manual
# PLEASE NOTE : simulation duration IS SET TO 0 ON PURPOSE - the appropriate duration is calculated below
pcomb1=[nedge,tres,sim_durn,cornerpos,iid_on_axis,mov_angle,flight_speed,call_durn,ipi_combs,hear_thresh,source_level,alpha,speedsound]
#pcomb2=[nedge,tres,sim_durn,cornerpos,iid_on_axis,mov_angle,flight_speed,call_durn,ipi_combs[1],hear_thresh,source_level,alpha,speedsound]
#pcomb3=[nedge,tres,sim_durn,cornerpos,iid_on_axis,mov_angle,flight_speed,call_durn,ipi_combs[2],hear_thresh,source_level,alpha,speedsound]
#pcomb4=[nedge,tres,sim_durn,cornerpos,iid_on_axis,mov_angle,flight_speed,call_durn,ipi_combs[3],hear_thresh,source_level,alpha,speedsound]
#pcomb5=[nedge,tres,sim_durn,cornerpos,iid_on_axis,mov_angle,flight_speed,call_durn,ipi_combs[4],hear_thresh,source_level,alpha,speedsound]
#

Nrep=1 # number of replicates to be run per parameter combination    
numcycles=4 # number of call cycles that will simulated for all parameter combinations


#function which calculates required simulation duration based on the number of call cycles to be studied:
#inputs:
# pulse interval (float), call duration (float), number of call cycles (integer),time resolution (float)
# output: number of iterations the simulation should be run for (integer)

def simdurcalc(PI,calldurn,numcallcycles,timeresn):
    simdurn=(numcallcycles+1)*(PI+calldurn) # duration of the simulation in seconds - the numcallcycles +1 is to include the bat that calls at the last momen before PI is up
    simdurn_iters=np.rint(simdurn/timeresn) # duration in integers round off number of iterations to closest integer
    return (simdurn_iters)
    


# combine the parameter combination in a list with sublists
Param_set=[pcomb1]

#ensure the correct simulation duration is calculated for each parameter combination:
for pcb in Param_set:
    pcb[2]=simdurcalc(pcb[8],pcb[7],numcycles,pcb[1]) # reassigning the value for number of iterations to be run
    

## status messages 
print ('simulation has started - and is running..... \n')
print (' \n %s Replicates per parameter combination with %s parameter combinations' %(Nrep,len(Param_set)) )
   


# begin the replicate simulations for each parameter combinations:

for repno in range(Nrep):
    print('\n.......the %s th replicate of all parameter combinations is now being run' % str(repno+1))
    for parcomb in Param_set:
            tgtdir=os.path.join(Results_dir,'Res%s'%str(repno)) # create a directory for each replicate of all pcombs run          
            onerun(tgtdir,parcomb)  # run the simulation for the pcomb               
                
   
print('\n all parameter combinations have been succesfully run - please go to %s to access the data. \n' %Results_dir)   
print (time.clock()-start) #uncomment to see how long the whole script takes 
    