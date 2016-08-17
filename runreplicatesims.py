# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 12:57:40 2016

@author: tbeleyur
"""
import sys
import os
import numpy as np
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
Results_dir='D:\\Bat_Project\\Res'


# import the script ('oneinst.py') which runs one single instance of the simulation : 
modulelocn='C:\\Users\\Claire\\Documents\\GitHub\\Bat-Project'

sys.path.append ( modulelocn ) # add the location of the module to the search path of python 
from oneinst import onerun  # import the one function from the module 


Nrep=100 # number of replicates to be run per parameter combination

# create the parameter combinations wanted : - right now it is purposefully manual
# PLEASE NOTE : simulation duration IS SET TO 0 ON PURPOSE - the appropriate duration is calculated below
pcomb1=[5,0.001,0,[1,1],2,0,5,0.003,0.080,-10,120,-1.7,340]
pcomb2=[5,0.001,0,[1,1],2,0,5,0.003,0.050,-10,120,-1.7,340]
pcomb3=[5,0.001,0,[1,1],2,0,5,0.003,0.040,-10,120,-1.7,340]
pcomb4=[5,0.001,0,[1,1],2,0,5,0.003,0.020,-10,120,-1.7,340]
pcomb5=[5,0.001,0,[1,1],2,0,5,0.003,0.010,-10,120,-1.7,340]




#function which calculates required simulation duration based on the number of call cycles to be studied:
#inputs:
# pulse interval (float), call duration (float), number of call cycles (integer),time resolution (float)
# output: number of iterations the simulation should be run for (integer)

def simdurcalc(PI,calldurn,numcallcycles,timeresn):
    simdurn=(numcallcycles+1)*(PI+calldurn) # duration of the simulation in seconds - the numcallcycles +1 is to include the bat that calls at the last momen before PI is up
    simdurn_iters=np.rint(simdurn/timeresn) # duration in integers round off number of iterations to closest integer
    return (simdurn_iters)
    
    
numcycles=9 # number of call cycles that will simulated for all parameter combinations


# combine the parameter combination in a list with sublists
Param_set=[pcomb1,pcomb2]

#ensure the correct simulation duration is calculated for each parameter combination:
for pcb in Param_set:
    pcb[2]=simdurcalc(pcb[8],pcb[7],numcycles,pcb[1]) # reassigning the value for number of iterations to be run
    



# begin the replicate simulations for each parameter combinations:

for repno in range(Nrep):
    for parcomb in Param_set:
            tgtdir=os.path.join(Results_dir,'Res%s'%str(repno)) # create a directory for each replicate of all pcombs run          
            onerun(tgtdir,parcomb)  # run the simulation for the pcomb               
                
   
   
print (time.clock()-start) #uncomment to see how long the whole script takes 
    