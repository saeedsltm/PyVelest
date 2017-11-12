import numpy as np
from datetime import datetime as dt
import pylab as plt
from os import rename, path
from matplotlib import rcParams
import sys
sys.path.append(path.join('tools'))
from PyNordicRW import Read_Nordic



plt.rc('font',family='Times New Roman')
plt.style.use('seaborn-paper')

rcParams['axes.labelsize'  ] = 15
rcParams['xtick.labelsize' ] = 14
rcParams['ytick.labelsize' ] = 14
rcParams['legend.fontsize' ] = 12
rcParams['legend.numpoints'] = 1

              

def PRD(inp_file='hyp.out',sta_inp='station.sta', falg_plot=False):

    '''
    A PRD function (PHASE-RESIDUAL-DISTANT) is weithing function based on
    three important parameters: 1- number of recorded phases [uses only P phase]
                                2- residuals
                                3- event-station distance
                                
    wich calculate the most important station when locating a cluster
    It uses three exponensial function in which each of them is related
    to one paramter(P,R,D) as describe bellow:

    F(p,r,d) = Exp(-[(pm-pi)/pm]) * Exp(-[(ri-rm)/rm]) * Exp(-[di/dm])

    pm = mean of all recorded phases for all stations
    pi = mean of all recorded phases for station (i)

    rm = mean of all residuals for all stations
    ri = mean of all residuals for station (i)

    dm = mean of all event-station distances for all stations
    di = mean of all event-station distance for station (i)

    the final weigth will be given to station with top value of F(p,r,d)
    
    '''

#___________________ Read a nordic file and calculate the reference station useing PRD weithing function [decribed below]

    s        = Read_Nordic(inp_file)
    sta      = {}
    used_sta = []
    
    with open(sta_inp) as f:

        for l in f:

            if '(' not in l:
                used_sta.append(l[:4].strip())

#__________ Prepare a dict as input for PRD function

    for i in s.keys():

        for j in s[i]['PHASE'].keys():

            if j not in sta.keys() and j in used_sta:

                sta[j] = {'pha':0,'res':[],'dis':[]}

            for phase in s[i]['PHASE'][j]['P']:
                
                if j in used_sta and s[i]['PHASE'][j]['P'][phase].has_key('TRES') and s[i]['PHASE'][j]['P'][phase].has_key('DIS') :

                    if s[i]['PHASE'][j]['P'][phase]['TRES'] != None and s[i]['PHASE'][j]['P'][phase]['DIS'] != None:

                        sta[j]['pha']+=1
                        sta[j]['res'].append(s[i]['PHASE'][j]['P'][phase]['TRES'])
                        sta[j]['dis'].append(s[i]['PHASE'][j]['P'][phase]['DIS'])
                    

    for p in sta.keys():

        if sta[p]['pha'] == 0:

            sta.pop(p)

    dic = sta
    
    mean_pha = []
    mean_res = []
    mean_dis = []
    non_val  = '     '

    for i in sorted(dic.keys()):

        for _ in range(dic[i]['res'].count(non_val)): dic[i]['res'].remove(non_val)
        for _ in range(dic[i]['dis'].count(non_val)): dic[i]['dis'].remove(non_val)

        mean_pha.append(dic[i]['pha'])
        mean_res+=dic[i]['res']
        mean_dis+=dic[i]['dis']

    mean_res = [i for i in mean_res if i]
    mean_dis = [i for i in mean_dis if i]

    mean_pha = np.mean(mean_pha)
    mean_res = 0.0              # Best value is zero :)
    mean_dis = np.mean(mean_dis)

    print "Mean Number of Phases: %5d # "%(mean_pha)
    print "Mean Distance        : %5d km"%(mean_dis)

    res = {}

    print "\nSTA  ** NumPhase ** MRes  ** MDis **  W\n"

    for i in sorted(dic.keys()):

        sta_res_mean = np.mean([k for k in dic[i]['res'] if k != None])
        sta_dis_mean = np.mean([k for k in dic[i]['dis'] if k != None])

        if not np.isnan(sta_res_mean.size) and not np.isnan(sta_dis_mean):

            f1 = np.exp(-(mean_pha - dic[i]['pha'])/mean_pha)
            f2 = np.exp(-abs(sta_res_mean))
            f3 = np.exp(-sta_dis_mean/mean_dis)
            w  = f1 * f2 * f3

            res[i] = w

            print "%-4s ** %5d    ** %6.2f ** %3d  ** %4.2f"%(i,dic[i]['pha'],sta_res_mean,sta_dis_mean,w)


    print '\n Refrence station =',max(res,key=res.get)
    print ''
    
    res = sorted(res.items(), key=lambda x: x[1], reverse=True)

    #_________SAVE TO station.sta
    
    inp_sta = {}
    tmp_sta = file('tmp_sta','w')

    with open('station.sta') as f:

        for l in f:

            if '(' in l:
            
                hed = l.strip()
                tmp_sta.write(hed+'\n')

            elif '(' not in l and l.strip():

                inp_sta[l[0:4].replace(' ','')]={'lat':l[4:12],'lon':l[14:22],'dep':l[23:27]}
    cc = len(res)
    for i in res:

        try:

            sta = i[0]

            tmp_sta.write("%-4s%7s  %8s %4s %-3s%2d  %1s   %-2s   %-2s                              \n"%(sta, inp_sta[sta]['lat'], inp_sta[sta]['lon'], inp_sta[sta]['dep'], "1", cc,"0.00", "0.00", "1"))

            cc-=1

        except KeyError:

            print '%s not in station.dat!'%(sta)
            
    tmp_sta.write('\n')
    tmp_sta.close()

    rename('tmp_sta','station_ref.sta')

    #_________PLOT

    if falg_plot == True:

        plt.figure(figsize=(16,8))
        plt.grid()

        for i in sorted(dic.keys()):

            if i == res[0][0]:

                plt.plot(dic[i]['dis'],dic[i]['res'],'o',label=i)
                
            else:

                plt.plot(dic[i]['dis'],dic[i]['res'],'o',label=i,alpha=.7)

        legend = plt.legend(numpoints=1,loc=4, ncol=4,fontsize=10, title='Station',shadow=True)
        frame = legend.get_frame()
        frame.set_facecolor(np.array([float(248)/float(255)]*3))
        plt.xlabel('Distance [km]', fontsize=18)
        plt.ylabel('Residuals [sec]', fontsize=18)
        plt.ylim(-2.0,2.0)
        plt.tight_layout() 
        plt.savefig('Ref_Station.png',dpi=400, bbox_inches='tight')

        print '\nSaving file "Ref_Station.png" ...\n'
        
#___________________Example run

#PRD(inp_file='hyp.out',sta_inp='station.sta', falg_plot=False)

import warnings
warnings.filterwarnings("ignore")
