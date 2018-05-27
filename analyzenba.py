# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 16:16:06 2017

@author: derick
"""

from __future__ import division, print_function
import numpy as np
import matplotlib as plt
import datetime
import time
import cPickle
import os
import copy


#==============================================================================

def readPrevSeasons(pathdfol, pathdump):
    
    seasons = os.listdir(pathdfol)
    # Remove other files, keep data folders only.
    # Arrange season folders, 1900s before 2000s
    seasons = (sorted([s for s in seasons if (int(s.split('-')[0])>50)]) + 
               sorted([s for s in seasons if (int(s.split('-')[0])<50)])   )
    seasonsched = {}
    teams = {}
    
    for regsea in seasons:
        fpath = pathdfol + '\\' + regsea
        monthsched = os.listdir(fpath)
        
        # Get game schedules for different months
        schedule = {}
        teams[regsea] = []
        
        for mo in monthsched:
            
            sched = np.loadtxt(fpath + '\\' + mo, dtype=str, delimiter=',',
                              usecols=range(6))
            sched = sched[1:]
            
            # Format schedule based on date and time
            for game in sched:
                date = datetime.date(*time.strptime(game[0],
                                                    "%a %b %d %Y")[:3])
                
                if date not in schedule.keys():
                    schedule[date] = []
                
                # Data before 2000-2001 season have different format                    
                if int(regsea.split('-')[0])<50:
                    
                    schedule[date].append(
                        [str(game[4]), -1 if not game[5] else int(game[5]),
                         str(game[2]), -1 if not game[3] else int(game[3])])
                    # Get teams
                    if game[4] not in teams[regsea]: 
                        teams[regsea].append(game[4])
                else:
                    schedule[date].append(
                        [str(game[3]), -1 if not game[4] else int(game[4]),
                         str(game[1]), -1 if not game[2] else int(game[2])])
                    # Get teams
                    if game[3] not in teams[regsea]: 
                        teams[regsea].append(game[3])
                
                
        seasonsched[regsea] = schedule
        
        print(regsea, 'Regular season done!')
    
    for filedump in ['seasons', 'seasonsched', 'teams']:
        dfile = open(pathdump + '\\' + filedump + '.p', 'wb')
        cPickle.dump(eval(filedump), dfile, protocol=cPickle.HIGHEST_PROTOCOL)
        dfile.close()
    
    
    print('\nDumped binary files:\nseasons.p\nseasonsched.p\nteams.p\n')
    
    return seasons, seasonsched, teams


def overview(schedseason):
    ''' dates as keys only, one season only '''
    
    # Get number of game dates
    dates = np.array(sorted(schedseason.keys()))
    ngames = sum([len(schedseason[gdate]) for gdate in schedseason])  
    
    # Identify dates with no games (these include the All-star weekend)
    sepdates = np.array([sep.days for sep in dates[1:] - dates[:-1]])
    noglen = sepdates[np.where(sepdates>1)] - 1
    nogames = dates[np.where(sepdates>1)] + datetime.timedelta(1)
    nogames = [[gdate + datetime.timedelta(n) for n in xrange(noglen[ind])] 
                                    for ind,gdate in enumerate(nogames)]
    
    if len(nogames)!=0:
        allstar =  [[i[0],i[-1]] for i in nogames 
                                        if i[0].month==2][0]
    else: allstar = [None,None]
    
    return dates, ngames, nogames, allstar

    
def summPrevSeasons(seasons, seasonsched, teams, pathdump, fname='schedsum'):

    report = ''
    allstar = {}
    
    for regsea in seasons:
        dates, ngames, nogames, allstar[regsea] = overview(seasonsched[regsea])
                        
        report += '{0} Regular Season\n'.format(regsea)
        report += 'Start of season: {0}\nEnd of season: {1}\n'.format(
                        dates[0].strftime('%x'), dates[-1].strftime('%x'))
        report += 'Number of game dates: {0}\n'.format(len(dates))
        report += 'Number of games: {0}\n'.format(ngames)
        report += 'Number of teams: {0}\n'.format(len(teams[regsea]))
        report += 'Games per team: {0}\n'.format(
                                            int(ngames*2/len(teams[regsea])))
        report += 'Dates without games: {}\n'.format(', '.join(
                                    ngame[0].strftime('%x') if len(ngame)==1
                                            else (ngame[0].strftime('%x') + 
                                            ' to ' + ngame[-1].strftime('%x'))
                                                        for ngame in nogames))
        
        if seasons.index(regsea) != len(seasons)-1: report += '\n'
    
    
    txtf = open(fname + '.txt', 'w')
    txtf.write(report)
    txtf.close()
    
    dfile = open(pathdump + '\\allstar.p', 'wb')
    cPickle.dump(allstar, dfile, protocol=cPickle.HIGHEST_PROTOCOL)
    dfile.close()    
    
    print('Schedule summary created: {0}.txt'.format(fname))
    
    return allstar

def datekey2teamkey(schedseason, teamlist):
    '''schedule for each team, input sched dates as keys'''

    # Get the schedule for each team
    teamsch = {t:[] for t in teamlist}
    
    for gdate in schedseason:
        for game in schedseason[gdate]:
            teamsch[game[0]].append([gdate,game,'Home'])
            teamsch[game[2]].append([gdate,game,'Road'])
    # Arrange the game datesv
    for t in teamsch:
        teamsch[t].sort()
    
    return teamsch

def teamkey2datekey(teamsch):
    ''' convert teams as keys to game dates as keys'''    
    
    schedseason = {}
    
    for t in teamsch:
        for game in teamsch[t]:
            if game[0] not in schedseason.keys():
                schedseason[game[0]] = []
            
            if game[1] not in schedseason[game[0]]:
                schedseason[game[0]].append(game[1])
    
    return schedseason
    

def teamsched(seasons, seasonsched, teams, pathdump):
    ''' Get schedule for teams.''' 

    schedteams = {}
    for regsea in seasons:
        
        schedteams[regsea] = datekey2teamkey(seasonsched[regsea],
                                             teams[regsea])
    
    dfile = open(pathdump + '\\schedteams.p', 'wb')
    cPickle.dump(schedteams, dfile, protocol=cPickle.HIGHEST_PROTOCOL)
    dfile.close()
    
    print('\nDumped binary file:\nschedteams.p\n')
    
    return schedteams

#=============================================================================




#==============================================================================

def cluster(arr,homeroad=False):
    count = 0
    label = np.zeros(len(arr),dtype=int)
    for i in xrange(len(arr)):
        if i == 0 and arr[i] == 1:
            count += 1
            label[i] = count
        elif arr[i]==1 and arr[i-1] == 1:
            label[i] = label[i-1]
        elif arr[i] == 1 and arr[i-1] == 0:
            count += 1
            label[i] = count
    return label

def b2b(label):
    '''edit to accommodate multiple back to backs to back ... UPDATE: ok na'''
    count = max(label)
    cplabel = np.copy(label)
    b2b = {}
    b2b2bloc = []
    b2b2balert = 0
    for c in xrange(1,count+1):
        if (label == c).sum() not in b2b: b2b[(label == c).sum()] = 1
        else:  b2b[(label == c).sum()] += 1
########## For locating (i.e. date) higher order back to backs in the season
#        if (label==c).sum()==2:
#            b2b2balert = 1
#            b2b2bloc.append(np.where(label==c))
#            cplabel -= c*(label==c)
##    if len(b2b2bloc)!=0:print b2b2bloc
#    if b2b2balert: ret = (b2b,b2b2bloc,cplabel)
#    else: ret = (b2b,)
#####################################################
    ret = b2b
    return ret
            
def winlose(team,data):
    teamhomrod = data.index(team)
    return 1 if data[teamhomrod+1]>data[teamhomrod-1] else 0

def diffday(diff):
    '''diff must be datetime.timedelta object'''
    return diff.days

diffday = np.vectorize(diffday)

def findB2B(gdates):
    '''gdates must be numpy array'''
    
    
    diffgdate = diffday(gdates[1:] - gdates[:-1])
    gamediff = np.concatenate((diffgdate==1,[False]))*1 
    #^ add False para pareho size with date list, for masking purposes,
        
    label = cluster(gamediff)
    b2bdetails = b2b(label)
    b2bdetails[0] = sum(diffgdate==0)   #0 gdate difference -> >1 game per day
    
    return b2bdetails

def findB2B_new(gdates):
    '''gdates must be numpy array
    gdates are integers from first day of seasons with day 0 as first day
    '''
    
    diffgdate = gdates[1:] - gdates[:-1]
    gamediff = np.concatenate((diffgdate==1,[False]))*1
    #^ add False para pareho size with date list, for masking purposes,
        
    label = cluster(gamediff)
    b2bdetails = b2b(label)
    b2bdetails[0] = sum(diffgdate==0)   #0 gdate difference -> >1 game per day
    
    return b2bdetails


def gcircdist(p1,p2):
    if p1==p2: return np.float64(0)
    p1 = np.float64(np.deg2rad(p1))
    p2 = np.float64(np.deg2rad(p2))
    centang = np.arccos(np.sin(p1[0])*np.sin(p2[0]) + 
                    np.cos(p1[0])*np.cos(p2[0])*np.cos(np.abs(p2[1]-p1[1])))
    return 6371*centang # distance in km
    

def teamDistance(tsched, coords):
    '''tsched schedule of team, list of team names played at the tnames home,
    coords with teams as keys, coords Location object (for now)
    
    tsched format: team starts an ends at home arena'''
    dist = 0
    for i in xrange(len(tsched)-1):
        if tsched[i] != tsched[i+1]:
            dist += gcircdist(coords[tsched[i]][1], coords[tsched[i+1]][1])
    return dist


#==============================================================================

def seasonDistance(seasched, coords):
    ''' input format: teamsched output 
    seasched teams as keys
    calculates all distances of all teams'''
    tdist = []
    for t in seasched:
        tsched = [t] + [g[1][0] for g in seasched[t]] + [t]
        tdist.append(teamDistance(tsched,coords))
    
    return sum(tdist),np.std(tdist)
    

    
def seasonDistanceV2(seasched, teams, coords, distlist, teamchanges):
    '''calculates only the distances for teams with changes
    list of teams with changes in schedule
    teams is a python list object
    DISTLIST must be ARRANGED according to `teams`
    '''
    distlistV2 = copy.copy(distlist)
    for t in teamchanges:
        tsched = [t] + [g[1][0] for g in seasched[t]] + [t]
        
        distlistV2[teams.index(t)] = teamDistance(tsched,coords)
        
    return distlistV2

def getTeams(gamelist):
    teamlist = []
    for game in gamelist:
        if game[0] not in teamlist: teamlist.append(game[0])
        if game[2] not in teamlist: teamlist.append(game[2])
        
    return teamlist

def cost1(seasched, coords, teams):
    seaschedteams = datekey2teamkey(seasched,teams)
    b2b = {}
    for t in seaschedteams:
        tb2b = findB2B(np.array([g[0] for g in seaschedteams[t]]))
        for b2border in tb2b:
            if b2border not in b2b.keys():
                b2b[b2border] = 0
            
            b2b[b2border] += tb2b[b2border]
    print(b2b)
    dist,_ = seasonDistance(seaschedteams,coords)

    b2bcost = sum([b2b[b2border] if b2border==1 
                        else 10000*float(b2b[b2border]) for b2border in b2b])
    
    cost = dist/100000. + b2bcost
    
    return cost

def cost1_new(seasched, coords, teams):
    seaschedteams = datekey2teamkey(seasched,teams)
    b2b = {}
    for t in seaschedteams:
        tb2b = findB2B_new(np.array([g[0] for g in seaschedteams[t]]))
        for b2border in tb2b:
            if b2border not in b2b.keys():
                b2b[b2border] = 0
            
            b2b[b2border] += tb2b[b2border]
#    print(b2b)
    dist,_ = seasonDistance(seaschedteams,coords)

    b2bcost = sum([b2b[b2border] if b2border==1 
                        else 10000*float(b2b[b2border]) for b2border in b2b])
    
    cost = dist/100000. + b2bcost
    
    return cost

################ For possible optimization of code for calculating B2B
################ i.e. calculating only for teamms with changes in schedule
##def getStats(seasched,coords,teams):
##    '''date keys are integers'''
##    
##    ## total dist per team list arranged according to `teams`
##

def makedistlist(seaschedteams, coords, teams):
    '''essentially the same as seasonDistance'''
    distlist = []
    for t in teams:
        tsched = [t] + [g[1][0] for g in seaschedteams[t]] + [t]
        distlist.append(teamDistance(tsched,coords))
    return distlist
    

def cost_optimized(seasched, coords, teams, distlist, b2b, b2blist,
                   glist=None, prevcost=0):
    '''generalized cost function optimized for performance
    
    input 
    seassched -> season schedle; integer key dict
    coords -> team name key dict (for now)
    teams -> team list
    
    ##############comp=True -> compare with previous cost -> true by default
    glist -> list of games with changes, compute cost only for these teams
    prevcost -> prevcost
    '''

    seaschedteams = datekey2teamkey(seasched,teams)
    
    tlist = getTeams(glist)
    
    distlistV2 = seasonDistanceV2(seaschedteams, 
                                  teams, coords, distlist, tlist)
    b2blistV2 = copy.copy(b2blist)
    b2bV2 = copy.copy(b2b)
    for t  in tlist:
        tb2b = findB2B_new(np.array([g[0] for g in seaschedteams[t]]))
        b2blistV2[teams.index(t)] = tb2b[1]
        b2bV2[t] = tb2b
    
    totalb2b={}
    for t in b2bV2:
        for b2border in b2bV2[t]:
            if b2border not in totalb2b.keys():
                totalb2b[b2border] = 0
        
            totalb2b[b2border] += b2bV2[t][b2border]
            

    b2bcost = sum([totalb2b[b2border] if b2border==1 
                else 2**b2border*10000*float(totalb2b[b2border]) for b2border in totalb2b])
    b2bspread = np.std(b2blistV2)
    
    cost = np.sum(distlistV2)/10000. + b2bcost + 10*b2bspread + np.std(distlistV2)/1000
    
    return cost, distlistV2, b2blistV2, b2bV2
        
def initCost(distlist,b2blist,b2b):

    totalb2b={}
    for t in b2b:
        for b2border in b2b[t]:
            if b2border not in totalb2b.keys():
                totalb2b[b2border] = 0
        
            totalb2b[b2border] += b2b[t][b2border]
            

    b2bcost = sum([totalb2b[b2border] if b2border==1 
                else 2**b2border*10000*float(totalb2b[b2border]) for b2border in totalb2b])
    b2bspread = np.std(b2blist)
    
    cost = np.sum(distlist)/10000. + b2bcost + 10*b2bspread + np.std(distlist)/1000
    
    return cost

def countViolations(seasched, teams, b2b, glist):
    '''count constraint violations, e.g. b2b2b, multiple games in a day
    
    returns totViolations,b2bV2, totalb2b
    '''

    seaschedteams = datekey2teamkey(seasched, teams)
    
    tlist = getTeams(glist)
    
    b2bV2 = copy.copy(b2b)
    for t  in tlist:
        tb2b = findB2B_new(np.array([g[0] for g in seaschedteams[t]]))
        b2bV2[t] = tb2b
    
    totalb2b={}
    for t in b2bV2:
        for b2border in b2bV2[t]:
            if b2border not in totalb2b.keys():
                totalb2b[b2border] = 0
        
            totalb2b[b2border] += b2bV2[t][b2border]
    
    Violations = np.array([ [ 2**b2border*100*totalb2b[b2border], totalb2b[b2border] ]
                          for b2border in totalb2b if b2border!=1 ])
    
    return Violations[:,0].sum(), Violations[:,1].sum(), b2bV2, totalb2b
    
    
def getCoordSeason(seas, currarenacoord, teamlist,arenaadjcoord):
    '''get list of coordinates for a given season `seas`
    assumes that adjustments list `arenaadjcoord` is loaded'''
    
    arenacoordadj = copy.deepcopy(currarenacoord)
    # starting year of season
    styear = (int(seas.split('-')[0]) + 2000 
                if int(seas.split('-')[0][0])<5 
                    else int(seas.split('-')[0]) + 1900)  
    # Check for adjustments
    for t in arenacoordadj.keys():
        # Check team changes
        if t not in teamlist:
            arenacoordadj.pop(t)
    # If there is adjustment for the year, change data
    for adj in arenaadjcoord:
        if styear in adj[0]:
            arenacoordadj[adj[1]] = adj[-1]   
    
    return arenacoordadj

#def saCost1_new(seasonsched):

#    plt.subplot(len(seasons),1,seasons.index(sea)+1)
#    plt.hist([len(seasonsched[sea][d]) for d in seasonsched[sea]])

#for sea in seasons:
#   print(np.std(([len(seasonsched[sea][d]) for d in seasonsched[sea]])))


    

