# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 12:57:40 2016

@author: tbeleyur
"""
import sys
import os
# script which runs multiple replicates of many parameter combinations:
# inputs : 
# results_dir: folder in which ALL the results will be stored
# Nrep : number of replicate simulations to be run per parameter set
# Param_set: list with sublists. [ [parcombi1],[parcombi2],parcombi3..] 

#outputs:
# results in the specificed folder - each parameter combination gets its own labelled folder within each 
# replicate numbered folder (Res#repnum)

#the target folder that you want all the results to be stored in : (folder doesn't need to exist ! ) 
Results_dir='C:\\Users\\tbeleyur\\Desktop\\final_results_TEST'


# import the script ('oneinst.py') which runs one single instance of the simulation : 
modulelocn='C:\\Users\\tbeleyur\\Google Drive\\Holger Goerlitz- IMPRS\\PHD_2015\\projects and analyses\\2016_jamming response modelling\\Coding\\Thejasvi codes and comments\\Bat-Project'

sys.path.append ( modulelocn ) # add the location of the module to the search path of python 
from oneinst import onerun  # import the one function from the module 


Nrep=5 # number of replicates to be run per parameter combination

# create the parameter combinations wanted : - right now it is purposefully manual
pcomb1=[3,0.001,30,[1,1],2,0,5,0.003,0.080,-10,120,-1.7]
pcomb2=[3,0.001,30,[1,1],2,0,5,0.007,0.080,-10,120,-1.7]


# combine the parameter combination in a list with sublists
Param_set=[pcomb1,pcomb2]


# begin the replicate simulations for each parameter combinations:

for repno in range(Nrep):
    for parcomb in Param_set:
            tgtdir=os.path.join(Results_dir,'Res%s'%str(repno)) # create a directory for each replicate of all pcombs run          
            onerun(tgtdir,parcomb)  # run the simulation for the pcomb               
                
    
    