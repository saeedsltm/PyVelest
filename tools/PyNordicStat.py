from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from numpy import zeros_like, ones_like
import os, sys
sys.path.append(os.path.join('utility'))
from PySTA0RW import Read_Nordic_Sta, Write_Nordic_Sta
from PyNordicRW import Read_Nordic
from initial_mpl import init_plotting_isi

"""
2017, Jun => Removed stations now displaying with a cross (x) in front of the station line in output.
"""

#___________________GET INPUT

#inp         = raw_input('\n\n+++ Nordic File?\n\n')
inp = 'select.out'
min_num_pha = raw_input('\n+++ Minimum number of recorded phases required in each station [def=10]:\n\n')
lon_rng     = raw_input('\n+++ Enter Longitude Min, Max (comma separated; def=Global):\n\n')
lat_rng     = raw_input('\n+++ Enter Latitude Min, Max (comma separated; def=Global):\n\n')
sta         = {}

if min_num_pha.split():

    min_num_pha = int(min_num_pha)

else:

    min_num_pha = 10
    
if lon_rng.split() and lat_rng.split():
    
    lon_rng = lon_rng.split(',')
    lat_rng = lat_rng.split(',')

    lon_min, lon_max = float(lon_rng[0]), float(lon_rng[1])
    lat_min, lat_max = float(lat_rng[0]), float(lat_rng[1])

else:

    lon_min, lon_max = -180, 180
    lat_min, lat_max =  -90,  90

#___________________READ NORDIC PHASE & STATION FILES

inp_nordic = Read_Nordic(inp)
sta0_orig  = Read_Nordic_Sta()

#___________________CHECK MINIMUM NUMBER OF PHASES REQUIRED IN EACH STATION.  

fi = file('statistic_report.out','w')
fi.write('Station Num_Pg Num_Pn Num_Sg Pol_C Pol_D   R_STA\n')
fi.write('*'*49+'\n')

Pg_sum    = 0
Pn_sum    = 0
Sg_sum    = 0
Pol_C_sum = 0
Pol_D_sum = 0
num_rm_s  = 0

nordic_sta = {}

for evt in inp_nordic.keys():

    for sta in inp_nordic[evt]['PHASE'].keys():

        if not nordic_sta.has_key(sta):

            nordic_sta[sta] = {'Num_Pg':0, 'Num_Pn':0, 'Num_Sg':0, 'Num_C':0, 'Num_D':0, }

        for p in inp_nordic[evt]['PHASE'][sta]['P']:

            if inp_nordic[evt]['PHASE'][sta]['P'][p]['POL'].upper() in ['C','+','U']:

                nordic_sta[sta]['Num_C']+=1

            if inp_nordic[evt]['PHASE'][sta]['P'][p]['POL'].upper() in ['D','-']:

                nordic_sta[sta]['Num_D']+=1
                
            if (p.upper() == 'PG' or p.upper() == 'P') and inp_nordic[evt]['PHASE'][sta]['P'][p]['WI'] != 4:

                nordic_sta[sta]['Num_Pg']+=1

            if p.upper() == 'PN' and inp_nordic[evt]['PHASE'][sta]['P'][p]['WI'] != 4:

                nordic_sta[sta]['Num_Pn']+=1

        for p in inp_nordic[evt]['PHASE'][sta]['S']:

            if (p.upper() == 'SG' or p.upper() == 'S') and inp_nordic[evt]['PHASE'][sta]['S'][p]['WI'] != 4:

                nordic_sta[sta]['Num_Sg']+=1

#___________________REMOVE UNUSED STATION FROM STATION0.HYP
                
for s in sta0_orig['STA'].keys():

    if not nordic_sta.has_key(s):

        del(sta0_orig['STA'][s])
        num_rm_s+=1

#___________________REMOVE STATIONS WHICH HAS NO MINIMUM QULIFICATION

for s in sorted(nordic_sta.keys()):

    num_pg  = nordic_sta[s]['Num_Pg']
    num_pn  = nordic_sta[s]['Num_Pn']
    num_sg  = nordic_sta[s]['Num_Sg']
    n_pol_d = nordic_sta[s]['Num_D']
    n_pol_c = nordic_sta[s]['Num_C']

    Pg_sum+=num_pg
    Pn_sum+=num_pn
    Sg_sum+=num_sg
    Pol_C_sum+=n_pol_c
    Pol_D_sum+=n_pol_d

    c1 = sta0_orig['STA'].has_key(s)
    
    if c1:
        
        c2 = (num_pg + num_pn + num_sg) < min_num_pha
        c3 = True
        c4 = True
        c3 = lon_min <= sta0_orig['STA'][s]['LON'] <= lon_max 
        c4 = lat_min <= sta0_orig['STA'][s]['LAT'] <= lat_max

        if (not c3 or not c4) or c2:

            del(sta0_orig['STA'][s])
            num_rm_s+=1
            fi.write('%-7s %6d %6d %6d %6d %6d   x\n'%(s, num_pg, num_pn, num_sg, n_pol_c, n_pol_d))

        else:

            fi.write('%-7s %6d %6d %6d %6d %6d   \n'%(s, num_pg, num_pn, num_sg, n_pol_c, n_pol_d))

fi.write('*'*42)
fi.write('\nTotal = %6d %6d %6d %6d %6d\n'%(Pg_sum, Pn_sum, Sg_sum, Pol_C_sum, Pol_D_sum))

#___________________WRITE OUTPUTS

fi.close()    
Write_Nordic_Sta(sta0_orig)

#___________________MAKE PLOT

init_plotting_isi(18,9)
plt.rcParams['axes.labelsize'] = 6
plt.rcParams['axes.titlesize'] = 6
    
ax  = plt.subplot(111, projection='3d')

x = []
y = []
z = []

for s in sta0_orig['STA'].keys():

    x.append(sta0_orig['STA'][s]['LON'])
    y.append(sta0_orig['STA'][s]['LAT'])
    z.append(nordic_sta[s]['Num_Pg']+nordic_sta[s]['Num_Pn']+nordic_sta[s]['Num_Sg'])

    ax.text(x[-1],y[-1],z[-1],s,ha='center',va='center')

xpos = x
ypos = y
zpos = zeros_like(x)
dx   = ones_like(x)*.02
dy   = ones_like(y)*.02
dz   = z

ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='red', zsort='average', alpha=.7)

ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_zlabel('Number of recorded phase (#)')

ax.xaxis._axinfo["grid"].update({"linewidth":.5, "color" : "gray", "linestyle":"--"})
ax.yaxis._axinfo["grid"].update({"linewidth":.5, "color" : "gray", "linestyle":"--"})
ax.zaxis._axinfo["grid"].update({"linewidth":.5, "color" : "gray", "linestyle":"--"})

ax.w_xaxis.set_pane_color((1.0, 1.0, 0.0, 0.5))
ax.w_yaxis.set_pane_color((1.0, 1.0, 0.0, 0.6))
ax.w_zaxis.set_pane_color((1.0, 1.0, 0.0, 0.7))

plt.tight_layout()
plt.savefig('statistical_plot.tiff')
plt.close()

#___________________MAKE SUMMARY

print ''
print '+++ Min No. of phases required in each station :',min_num_pha
print '+++ Min,Max of Longitude for selecting station :',lon_min, lon_max
print '+++ Min,Max of Latitude for selecting station  :',lat_min, lat_max
print '+++ Number of removed stations                 :',num_rm_s
print '+++ Statistical Output                         : statistic_report.out'
print '+++ New station Output                         : STATION0_NEW.HYP'
