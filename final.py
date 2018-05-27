# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 16:51:42 2017

@author: derick
"""

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import cPickle
import analyzenba as analnba
import os
import time
import datetime
import copy

# Get current working directory
pathmain = os.getcwd()

# Check for data
orgdatapath = pathmain + '\\Organized_Data'

#==============================================================================
# if ('Organized_Data' not in os.listdir(pathmain)): os.mkdir(orgdatapath)
# 
# if 'seasonsched.p' not in os.listdir(orgdatapath):
#     print('Previous schedules not yet read!\nOrganizing data...')
#     
#     seasons, seasonsched, teams = analnba.readPrevSeasons(pathmain + \
#                                                         '\\Raw_Data',
#                                                         orgdatapath)
#     
#     print('Data organized!\n')
#     
# else:
#     print('Data organized!\nExtracting data...')
#     
#     # Load seasons data
#     dfile = open(orgdatapath + '\\seasons.p', 'rb')
#     seasons = cPickle.load(dfile)
#     dfile.close()
#     
#     # Load season schedules data
#     dfile = open(orgdatapath + '\\seasonsched.p', 'rb')
#     seasonsched = cPickle.load(dfile)
#     dfile.close()
#     
#     # Load team name/participation data
#     dfile = open(orgdatapath + '\\teams.p', 'rb')
#     teams = cPickle.load(dfile)
#     dfile.close()
#     
#     print('Data extraction done!')
# 
# 
# # Summary of schedules
# if not ('schedsum.txt' in os.listdir(pathmain) and 
#         'allstar.p' in os.listdir(orgdatapath)):
#     print('Creating summary of previous schedules ...')
#     allstar = analnba.summPrevSeasons(seasons, seasonsched, teams, orgdatapath)
#     
# else:
#     print('Reading all-star data ...')
#     dfile = open(orgdatapath + r'\allstar.p', 'rb')
#     allstar = cPickle.load(dfile)
#     dfile.close()
# 
# if 'schedteams.p' not in os.listdir(orgdatapath):
#     print('Reading schedules per team...')
#     schedteams = analnba.teamsched(seasons, seasonsched, teams, orgdatapath)
# else:
#     # Load schedule of teams data
#     dfile = open(orgdatapath + r'\schedteams.p', 'rb')
#     schedteams = cPickle.load(dfile)
#     dfile.close()
#     
#     print('Schedule of teams loaded!')
# 
# 
# dfile = open(orgdatapath + r'\arenacoord1516.p', 'rb')
# arenacoordnow = cPickle.load(dfile)
# dfile.close()
# 
# dfile = open(orgdatapath + r'\arenaadjust.p', 'rb')
# arenaadjcoord = cPickle.load(dfile)
# dfile.close()
#==============================================================================

#==============================================================================
# ls = ['cgdANDgds','cgdANDbds','all','gdsANDbds']
# maxiter = 300000
# 
# pathdump = pathmain + '\\locSearchAnal3'
# report = ''
# 
# for searchmeth in ls:
#     # Do 10 times
#     tstart = time.time()
#     locseadata = []
#     
#     for itercount in xrange(10):
#         print(searchmeth,itercount)
# 
#         # Removes violations via random search greedy; test local search moves
#         #tstart = time.time()
#         
#         # Initialize schedule; random initialization
#         regsea = seasons[-1] # 17-18 regular season
#         gamepool = [g for gdate in seasonsched[regsea] for g in seasonsched[regsea][gdate]]
#         gamepool = map(list, np.random.permutation(gamepool))
#         presched = {}
#         count = 0
#         
#         gamedates = np.array(sorted(seasonsched[regsea].keys()))
#         gdatemark = np.cumsum(analnba.diffday(gamedates[1:] - gamedates[:-1]))
#         gdatemark = np.concatenate(([0],gdatemark))
#         
#         for gdate in gdatemark:
#             n = np.random.randint(6,9)
#             presched[gdate] = gamepool[count:count+n]
#             count += n
#         for n in xrange(1230-count):       ##### TO OPTIMIZE
#             presched[np.random.choice(gdatemark)].append(gamepool[count+n])
#         
#         preschedteams = analnba.datekey2teamkey(presched,teams[regsea])
#         
#         
#         b2b = {}
#         totalb2b = {}
#         for t in teams[regsea]:
#             tb2b = analnba.findB2B_new(np.array([g[0] for g in preschedteams[t]]))
#             b2b[t] = tb2b
#             
#             for b2border in tb2b:
#                 if b2border not in totalb2b.keys():
#                     totalb2b[b2border] = 0
#                 
#                 totalb2b[b2border] += tb2b[b2border]
#         
#         
#         Violations = np.array([ [ 2**b2border*100*totalb2b[b2border], totalb2b[b2border] ]
#                                   for b2border in totalb2b if b2border!=1 ])
#         costViolations = Violations[:,0].sum()
#         numViolations = Violations[:,1].sum()
#         
#         #        T = 100.
#         count = 0
# #        print('\t', costViolations, numViolations)
#         
#         consviolate = []
#         
#         consviolate.append([costViolations, numViolations])
#         
#         
#         while costViolations>0 and count<maxiter:
#             if count%100 == 0:
#         #                T = 0.999*T
#                 print('\t', costViolations, numViolations)
#         #    cost = analnba.cost1_new(presched,arenacoordnow,teams[regsea])
#             
#             mod = np.random.choice(presched.keys(),2,replace=0)
#             
#                  
#             
#             if searchmeth=='cgd': search=3
#             elif searchmeth=='cgdANDgds': search = np.random.randint(2,4)
#             elif searchmeth=='cgdANDbds': search = np.random.choice([1,3])
#             elif searchmeth=='all': search = np.random.randint(1,4)
#             elif searchmeth=='gdsANDbds': search = np.random.randint(1,3)
#             elif searchmeth=='gds': search=2
#             else : search=1
#             
#                 
#             if search==1:
#                 presched[mod[0]], presched[mod[1]] = presched[mod[1]], presched[mod[0]]
#                 
#                 glist = presched[mod[0]] + presched[mod[1]] 
#         
#             if search==2:
#                 sw0 = np.random.randint(len(presched[mod[0]]))
#                 sw1 = np.random.randint(len(presched[mod[1]]))
#                 
#                 presched[mod[0]][sw0], presched[mod[1]][sw1] = \
#                             presched[mod[1]][sw1], presched[mod[0]][sw0]
#                 
#                 glist = [ presched[mod[0]][sw0], presched [mod[1]][sw1] ]
#         #        temp = presched[mod[0]].pop(sw0)
#         #        presched[mod[1]].append(temp)
#         #        presched[mod[0]].append(presched[mod[1]].pop(sw1))
#             
#             if search==3:
#                 while len(presched[mod[0]]) <4:
#                     temp  = presched.keys()
#                     temp.pop(temp.index(mod[1]))
#                     mod[0] = np.random.choice(temp)
#                     
#                 sw0 = np.random.randint(len(presched[mod[0]]))
#                 
#                 temp = presched[mod[0]].pop(sw0)
#                 presched[mod[1]].append(temp)
#                 
#                 glist = [ temp ]
#             
#             
#             
#         #    costnew = analnba.cost1_new(presched,arenacoordnow,teams[regsea])
#         #    costnew, distlistV2, b2blistV2, b2bV2 = analnba.cost_optimized(presched, 
#         #                                    arenacoordnow, teams[regsea],
#         #                                    distlist, b2b, b2blist, glist)
#             costVio,numVio,b2bV2,totalb2bV2 = analnba.countViolations(presched,
#                                                               teams[regsea],b2b,glist)
#             
#             
#             costDiff = costVio - costViolations
#         #    print(costdiff)
#             
#             if costDiff>0:
#         #                if np.random.random()<np.exp(-costDiff/T): 
#         #                    costViolations = costVio
#         #                    numViolations = numVio
#         #                    b2b = b2bV2
#         #                else:            
#                 if search==1:
#                     presched[mod[0]], presched[mod[1]] = presched[mod[1]], presched[mod[0]]
#             
#                 if search==2:
#                     
#                     
#                     presched[mod[0]][sw0], presched[mod[1]][sw1] = \
#                                 presched[mod[1]][sw1], presched[mod[0]][sw0]
#                 
#                 if search==3:                
#                     temp = presched[mod[1]].pop(-1)
#                     presched[mod[0]].append(temp)
#                         
#             else:
#                 costViolations = costVio
#                 numViolations = numVio
#                 b2b = b2bV2
#             
#                 totalb2b = totalb2bV2
#             
#             count += 1
#             
#             consviolate.append([costViolations, numViolations])
#         
#         locseadata.append(consviolate)
#         
#     tend = time.time()
#     
#     dfile = open(pathdump + '\\%s.p' % searchmeth, 'wb')
#     cPickle.dump(locseadata, dfile, protocol=cPickle.HIGHEST_PROTOCOL)
#     dfile.close()
#     
#     report += '%-10s \t' % searchmeth
#     report += str((tend-tstart)/60.)
#     report += '\n'
#     
#     txtf = open(pathdump + '\\timeReport.txt', 'w')
#     txtf.write(report)
#     txtf.close()
#     
# #        print((tend-tstart)/60)
# 
# os.system('shutdown -s -t 180')
# 
#==============================================================================


# General annealing all variables and constraints at the same time
tstart = time.time()
#np.random.set_state(state)
regsea = '03-04' #seasons[-1] # 17-18 regular season

gamepool = [g for gdate in seasonsched[regsea] for g in seasonsched[regsea][gdate]]

gamepool = map(list,np.random.permutation(gamepool))

presched = {}
count = 0

gamedates = np.array(sorted(seasonsched[regsea].keys()))
gdatemark = np.cumsum(analnba.diffday(gamedates[1:] - gamedates[:-1]))
gdatemark = np.concatenate(([0],gdatemark))

for gdate in gdatemark:
    n = np.random.randint(6,9)
    presched[gdate] = gamepool[count:count+n]
    count += n
for n in xrange(1189-count):       ##### TO OPTIMIZE ## 1189 before 2004, 1230 after
    presched[np.random.choice(gdatemark)].append(gamepool[count+n])

preschedteams = analnba.datekey2teamkey(presched,teams[regsea])

#b2b = {}
#for t in preschedteams:
#    tb2b = analnba.findB2B_new(np.array([g[0] for g in preschedteams[t]]))
#    for b2border in tb2b:
#        if b2border not in b2b.keys():
#            b2b[b2border] = 0
#        
#        b2b[b2border] += tb2b[b2border]

b2blist = []
b2b = {}
for t in teams[regsea]:
    tb2b = analnba.findB2B_new(np.array([g[0] for g in preschedteams[t]]))
    b2blist.append(tb2b[1])
    b2b[t] = tb2b

totalb2b = {}
for t in b2b:
    for b2border in b2b[t]:
        if b2border not in totalb2b.keys():
            totalb2b[b2border] = 0
        
        totalb2b[b2border] += b2b[t][b2border]

arenacoordthen = analnba.getCoordSeason(regsea,arenacoordnow,teams[regsea],arenaadjcoord)

distlist = analnba.makedistlist(preschedteams,arenacoordthen,teams[regsea])

cost = analnba.initCost(distlist, b2blist, b2b)

#cost = analnba.cost1_new(presched,arenacoordnow,teams[regsea])
#coststart=cost
#
T = 1000.
count = 0
print(T,cost)
costlist=[]
distparams = []
b2bparams = []
consviolate = []
#costdiff=1e20
mincost = 1e50
#
costlist.append(cost)
distparams.append([np.sum(distlist), np.std(distlist)])
b2bparams.append([np.sum(b2blist), np.std(b2blist)])
consviolate.append(sum(totalb2b[b2border] for b2border in totalb2b if b2border != 1))
##rollstd = np.std(costlist[count+5-1-5:count+5-1])
##searchlist = np.random.randint(1,4,2000)
while count<1000000:#cost>1000000:
    if count%100 == 0 and count>0:
        T = 0.999*T
        print(T,cost)
#    cost = analnba.cost1_new(presched,arenacoordnow,teams[regsea])
    
    mod = np.random.choice(presched.keys(),2,replace=0)
    
    
#    if count < 10000:
#    search = np.random.randint(1,4)
#    else:
    search = np.random.randint(2,4)
        
    if search==1:
        presched[mod[0]], presched[mod[1]] = presched[mod[1]], presched[mod[0]]
        
        glist = presched[mod[0]] + presched[mod[1]]

    if search==2:
        sw0 = np.random.randint(len(presched[mod[0]]))
        sw1 = np.random.randint(len(presched[mod[1]]))
        
        presched[mod[0]][sw0], presched[mod[1]][sw1] = \
                    presched[mod[1]][sw1], presched[mod[0]][sw0]
        
        glist = [ presched[mod[0]][sw0], presched[mod[1]][sw1] ]
#        temp = presched[mod[0]].pop(sw0)
#        presched[mod[1]].append(temp)
#        presched[mod[0]].append(presched[mod[1]].pop(sw1))
    
    if search==3:
        while len(presched[mod[0]]) <4:
            temp  = presched.keys()
            temp.pop(temp.index(mod[1]))
            mod[0] = np.random.choice(temp)
            
        sw0 = np.random.randint(len(presched[mod[0]]))
        
        temp = presched[mod[0]].pop(sw0)
        presched[mod[1]].append(temp)
        
        glist = [ temp ]
    
    
    
#    costnew = analnba.cost1_new(presched,arenacoordnow,teams[regsea])
    costnew, distlistV2, b2blistV2, b2bV2 = analnba.cost_optimized(presched, 
                                    arenacoordthen, teams[regsea],
                                    distlist, b2b, b2blist, glist)
    
    if costnew<mincost:
        bestsched = copy.deepcopy(presched)
        mincost = costnew
    
    costdiff = costnew-cost
#    print(costdiff)
    
    if costnew>cost:
        if np.random.random()<np.exp(-costdiff/T): 
            cost = costnew
            distlist = distlistV2
            b2blist = b2blistV2
            b2b = b2bV2
        else:            
            if search==1:
                presched[mod[0]], presched[mod[1]] = presched[mod[1]], presched[mod[0]]
        
            if search==2:
                
                
                presched[mod[0]][sw0], presched[mod[1]][sw1] = \
                            presched[mod[1]][sw1], presched[mod[0]][sw0]
            
            if search==3:                
                temp = presched[mod[1]].pop(-1)
                presched[mod[0]].append(temp)
                
    else:
        cost=costnew
        distlist = distlistV2
        b2blist = b2blistV2
        b2b = b2bV2
    
    totalb2b = {}
    for t in b2b:
        for b2border in b2b[t]:
            if b2border not in totalb2b.keys():
                totalb2b[b2border] = 0
            
            totalb2b[b2border] += b2b[t][b2border]
    
    count += 1
    costlist.append(cost)
    distparams.append([np.sum(distlist), np.std(distlist)])
    b2bparams.append([np.sum(b2blist), np.std(b2blist)])
    consviolate.append(sum(totalb2b[b2border] for b2border in totalb2b if b2border != 1))
    
    tend = time.time()
#    rollstd = np.std(costlist[count+5-5-1:count+5-1])
#costlist.append(cost)
tend = time.time()
print((tend-tstart)/60)


