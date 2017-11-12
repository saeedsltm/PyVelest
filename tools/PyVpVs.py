import os
import sys
sys.path.append(os.path.join('utility'))
from PyNordicRW import Read_Nordic
from datetime import datetime as dt
import operator
from numpy import array, sqrt
from scipy.stats import linregress
import pylab as plt
from statsmodels.formula.api import ols
from initial_mpl import init_plotting_isi

phase_dic = {}

def outlier(x, y, r=0.5):

    regression = ols("data ~ x", data=dict(data=y, x=x)).fit()
    test       = regression.outlier_test()
    outliers_ind   = list((i for i,t in enumerate(test.iloc[:,2]) if t < r))
    good_id        = list(set(range(len(x)))-set(outliers_ind))

    return good_id, outliers_ind


def vpvs(d):

    for evt in sorted(d.keys()):
        
        a        = d[evt]
        ot_year  = a['HEADER']['L1']['Year']
        ot_month = a['HEADER']['L1']['Month']
        ot_day   = a['HEADER']['L1']['Day']
        ot_hour  = a['HEADER']['L1']['Hour']
        ot_min   = a['HEADER']['L1']['Min']
        ot_sec   = a['HEADER']['L1']['Sec']
        ot_msec  = int((ot_sec - int(ot_sec))*1e6)

        if ot_sec >= 60:

            ot_sec = ot_sec - 60
            ot_min+=1

        if ot_min >= 60:

            ot_min = ot_min - 60
            ot_hour+=1

        if ot_hour >= 24:

            ot_hour = ot_hour - 24
            ot_day+=1
                
        ot = dt(ot_year,ot_month,ot_day,ot_hour,ot_min,int(ot_sec),ot_msec)

        if not phase_dic.has_key(ot):

            phase_dic[ot] = {}
                    
        for i in a['PHASE']:

            if not phase_dic[ot].has_key(i):

                phase_dic[ot][i] = {'P':None, 'S':None}
                    
            otday    = ot_day

            for Pphase in a['PHASE'][i]['P']:
            
                try:

                    ar_hour = a['PHASE'][i]['P'][Pphase]['Hour']
                    ar_min  = a['PHASE'][i]['P'][Pphase]['Min']
                    ar_sec  = a['PHASE'][i]['P'][Pphase]['Sec']
                    ar_msec = int((ar_sec - int(ar_sec))*1e6)

                except (TypeError, ValueError):
                   
                    break

                if ar_sec >= 60:

                    ar_sec = ar_sec - 60
                    ar_min+=1

                if ar_min >= 60:

                    ar_min = ar_min - 60
                    ar_hour+=1

                if ar_hour >= 24:

                    ar_hour = ar_hour - 24
                    otday+=1
                        
                ar    = dt(ot_year,ot_month,otday,ar_hour,ar_min,int(ar_sec),ar_msec)  
                delta = ar - ot
                delta = delta.total_seconds()

                if delta > 0.0 and a['PHASE'][i]['P'][Pphase]['WI'] != 4:

                    phase_dic[ot][i]['P'] = delta
                    

            otday    = ot_day
            
            for Sphase in a['PHASE'][i]['S']:

                ar_hour = a['PHASE'][i]['S'][Sphase]['Hour']
                ar_min  = a['PHASE'][i]['S'][Sphase]['Min']
                ar_sec  = a['PHASE'][i]['S'][Sphase]['Sec']
                ar_msec = int((ar_sec - int(ar_sec))*1e6)
              
                if ar_sec >= 60:

                    ar_sec = ar_sec - 60
                    ar_min+=1

                if ar_min >= 60:

                    ar_min = ar_min - 60
                    ar_hour+=1

                if ar_hour >= 24:

                    ar_hour = ar_hour - 24
                    otday+=1

                ar    = dt(ot_year,ot_month,otday,ar_hour,ar_min,int(ar_sec),ar_msec)
                delta = ar - ot
                delta = delta.total_seconds()
                
                if delta > 0.0 and a['PHASE'][i]['S'][Sphase]['WI'] != 4:

                    phase_dic[ot][i]['S'] = delta

            if not phase_dic[ot][i]['S'] or not phase_dic[ot][i]['P']:

                del phase_dic[ot][i]


    # Pepare for plot
    
    p_time = []
    s_time = []

    for ev in sorted(phase_dic):

        if phase_dic[ev]:

            sorted_ev = sorted(phase_dic[ev].items(), key=operator.itemgetter(1))

            x = sorted_ev[0][1]['P']
            y = sorted_ev[0][1]['S']
                
            for p in sorted_ev:

                p_time.append(p[1]['P'] - x)
                s_time.append(p[1]['S'] - y)
    
    googd_id, outliers_ind = outlier(p_time,s_time)
    god_p        = []
    god_s        = []
    out_p        = []
    out_s        = []
    
    for _ in googd_id:

        god_p.append(p_time[_])
        god_s.append(s_time[_])

    for _ in outliers_ind:

        out_p.append(p_time[_])
        out_s.append(s_time[_])

    slope, intercept, r_value, p_value, std_err = linregress(god_p, god_s)

    y = slope*array(god_p) + intercept


    init_plotting_isi(9,6)

    ax = plt.subplot(111)

    ax.plot(god_p, god_s, marker='o', ls='', ms=3, mew=.1, mfc='r', mec='k')
    ax.plot(god_p, y, marker='', ls='--', lw=1, color='b')
    ax.scatter(out_p,out_s, marker='x', c='k', s=5, linewidth=.75, zorder=100)
    ax.set_xlim(0,max(p_time))
    ax.set_ylim(0,max(s_time))
    ax.set_xlabel('$P_{j}-P_{i}$ (s)')
    ax.set_ylabel('$S_{j}-S_{i}$ (s)')
    ax.set_title('Vp/Vs=%4.2f; N=%d; CC=%.2f; StdErr=%.2f'%(slope,len(god_p),r_value, std_err), fontsize=6)
    ax.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)

    plt.tight_layout()
    plt.savefig('VpVs_%4.2f.tiff'%(slope),dpi=300)
    plt.close()
    print '\n+++ Vp/Vs=%.2f\n'%(slope)


inp = raw_input('\n\n+++ Enter NORDIC file name:\n\n')

vpvs(Read_Nordic(inp))
