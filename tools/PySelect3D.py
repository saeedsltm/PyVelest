from numpy import genfromtxt, sqrt, arange, pi, meshgrid, linspace, finfo, array, mean
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import operator, sys, os
sys.path.append(os.path.join('utility'))
from PyNordicRW import Read_Nordic, Write_Nordic
from initial_mpl import init_plotting_isi
from matplotlib import ticker, gridspec




if not os.path.exists('select3D.par'):

    print '\n+++ Warning!\n'
    print '- No "select.par" file was found in current directory!'
    print '- Default values will be used and you can change them in the next run.\n'

    f = open('select3D.par','w')
    f.write("""
#
# Parameter file for selection program (PySelect),
# in a 3D space. Set the following parameters acc-
# ording to your study region. Avoide defining 
# small cells othervise it may take a long time!
#
#____________________________
#
LAT_MIN            = 35.40        # Minimum Latitude.
LAT_MAX            = 36.00        # Minimum Latitude.
LON_MIN            = 51.40        # Minimum Latitude.
LON_MAX            = 52.90        # Minimum Latitude.
DEP_MAX            = 40.00        # Maximum depth.
H_CELL_SIZE        = 20.00        # Cell size in horizontal [km^2].
Z_CELL_SIZE        =  5.00        # Cell size in depth [km].
NUM_EVT_MIN        = 1            # Minimum number of event in each cell.
NUM_EVT_MAX        = 15           # Maximum number of event in each cell. Type "A" for automatic detection or specefy a number.
Gap_Max            = 180 , 1.0    # Maximum Azimuthal Gap allowded for each event to be selected and it's weight.
NUM_STA_MIN        = 5   , 1.0    # Minumum number of station for each event to be selected and it's weight.
RMS_MAX            = 0.5 , 1.0    # Maximum RMS allowded for each event to be selected and it's weight.
HER_MAX            = 5.0 , 1.0    # Maximum horizontal error for each event to be selected and it's weight.
ZER_MAX            = 5.0 , 1.0    # Maximum depth error for each event to be selected and it's weight.
STA_MDIS           = 15.0, 1.0    # Station with minimum distance for event to be selected and it's weight.
EVENT_W            = 50           # Percentage of mean of all event's weight defines the Effective Weight. Events with weight more than Effective Weight, will be selected.
#
""")
    
    f.close()

select_par  = genfromtxt('select3D.par', delimiter='=', dtype=str)

#___________________SOME USEFULL FUNCTIONS

def k2d(kilometer, radius=6371):

    return kilometer / (2.0 * radius * pi / 360.0)

def d2k(degrees, radius=6371):

    return degrees * (2.0 * radius * pi / 360.0)

def check(v):

    if v>1:
        
        return 1.0
    
    elif v<0:
        
        return 0.0

    else:

        return v

def get_cells(sel_dic, d=2):

    X = []
    Y = []
    Z = []

    if d == 2:
        
        for c in sel_dic:

            x = sel_dic[c]['geo']['x']
            y = sel_dic[c]['geo']['y']

            x = [x[0],x[0],x[1],x[1],x[0]]
            y = [y[0],y[1],y[1],y[0],y[0]]

            X.append(x)
            Y.append(y)

        return X,Y

    if d == 3:

        for c in sel_dic:

            x = sel_dic[c]['geo']['x']
            y = sel_dic[c]['geo']['y']
            z = sel_dic[c]['geo']['z']

            z = list(-array(z))

            x = [x[0],x[0],x[0],x[1],x[1],x[1],x[1],x[0], x[0],x[0],
                 x[0],x[0],x[0],x[1],x[1],x[1],x[1],x[0], x[0],x[0],
                 x[0],x[0],x[0],x[1],x[1],x[1],x[1],x[1],x[1],x[1]]
            
            y = [y[0],y[0],y[0],y[0],y[0],y[0],y[0],y[0], y[0],y[1],
                 y[1],y[1],y[1],y[1],y[1],y[1],y[1],y[1], y[1],y[1],
                 y[1],y[0],y[0],y[0],y[0],y[1],y[1],y[1],y[1],y[0]]
            
            z = [z[0],z[1],z[0],z[0],z[0],z[1],z[1],z[1], z[1],z[1],
                 z[0],z[1],z[0],z[0],z[0],z[1],z[1],z[1], z[1],z[0],
                 z[0],z[0],z[0],z[0],z[0],z[0],z[0],z[1],z[1],z[1]]
            
            X.append(x)
            Y.append(y)
            Z.append(z)

        return X,Y,Z

#___________________READ SELECTION PARAMAETERS

inp     = raw_input('\n\n+++ Enter NORDIC file name:\n\n')

par_dic = {}

print '\n+++ Start selecting event(s) in a 3D space, Please Wait...\n'

for i in select_par:

    try:

        par_dic[i[0].strip()] = float(i[1])

    except ValueError:

        par_dic[i[0].strip()] = i[1].strip()

gap_max        = float(par_dic['Gap_Max'].split(',')[0])
gap_max_wt     = float(par_dic['Gap_Max'].split(',')[1])
num_sta_min    = float(par_dic['NUM_STA_MIN'].split(',')[0])
num_sta_min_wt = float(par_dic['NUM_STA_MIN'].split(',')[1])
rms_max        = float(par_dic['RMS_MAX'].split(',')[0])
rms_max_wt     = float(par_dic['RMS_MAX'].split(',')[1])
her_max        = float(par_dic['HER_MAX'].split(',')[0])
her_max_wt     = float(par_dic['HER_MAX'].split(',')[1])
zer_max        = float(par_dic['ZER_MAX'].split(',')[0])
zer_max_wt     = float(par_dic['ZER_MAX'].split(',')[1])
sta_mdis       = float(par_dic['STA_MDIS'].split(',')[0])
sta_mdis_wt    = float(par_dic['STA_MDIS'].split(',')[1])
evt_wt         = float(par_dic['EVENT_W'])
lat_min        = par_dic['LAT_MIN']
lat_max        = par_dic['LAT_MAX']
lat_len        = lat_max - lat_min
lon_min        = par_dic['LON_MIN']
lon_max        = par_dic['LON_MAX']
lon_len        = lon_max - lon_min
data           = Read_Nordic(inp)
out_dic        = {}
w_mean         = 0

#___________________GENERATE CELL GEOMETRIES

sta_list = []

for e in data:

    for s in data[e]['PHASE']:

        sta_list.append(s)

sta_list = sorted(set(sta_list))

sel_dic = {}

X = arange(par_dic['LON_MIN'], par_dic['LON_MAX'], k2d(sqrt(par_dic['H_CELL_SIZE'])))
Y = arange(par_dic['LAT_MIN'], par_dic['LAT_MAX'], k2d(sqrt(par_dic['H_CELL_SIZE'])))
Z = arange(0.0               , par_dic['DEP_MAX'], par_dic['Z_CELL_SIZE'])

cell_size = par_dic['H_CELL_SIZE'] * par_dic['Z_CELL_SIZE']
num_cell  = (len(X)-1)*(len(Y)-1)*(len(Z)-1)


if par_dic['NUM_EVT_MAX'] == 'A':

    max_nevt  = len(data)/num_cell

else:

    max_nevt = int(par_dic['NUM_EVT_MAX'])

c = 1

for z in range(len(Z)-1):

    for x in range(len(X)-1):

        for y in range(len(Y)-1):

            sel_dic[c] = {'geo':{'x':(X[x], X[x+1]),'y':(Y[y], Y[y+1]) ,'z':(Z[z], Z[z+1])},'evt_list':{}}
            c+=1
    
#___________________EXTRACT EVENT INTO CELLS

x_ini, x_fin = [], []
y_ini, y_fin = [], []
z_ini, z_fin = [], []
tot_evt_chk  = 0

for evt in data:

    try:
        
        x_ini.append(float(data[evt]['HEADER']['L1']['Lon']))
        y_ini.append(float(data[evt]['HEADER']['L1']['Lat']))
        z_ini.append(-float(data[evt]['HEADER']['L1']['Dep']))

        lat = data[evt]['HEADER']['L1']['Lat']
        lon = data[evt]['HEADER']['L1']['Lon']
        dep = data[evt]['HEADER']['L1']['Dep']
        gap = data[evt]['HEADER']['LE']['GAP']
        nst = data[evt]['HEADER']['L1']['NumS']
        rms = data[evt]['HEADER']['L1']['RMS']
        her = sqrt(data[evt]['HEADER']['LE']['LAE']**2 + data[evt]['HEADER']['LE']['LOE']**2) 
        zer = data[evt]['HEADER']['LE']['DEE']
        smd = []
    
    except (TypeError, ValueError):

        #print '+++ Skipping Event! Check Header =>', evt
        tot_evt_chk+=1

        continue

    for _ in data[evt]['PHASE']:

        for __ in data[evt]['PHASE'][_]['P']:

            if 'DIS' in data[evt]['PHASE'][_]['P'][__]:

                smd.append(data[evt]['PHASE'][_]['P'][__]['DIS'])

    try:
        smd = min(smd)
    except:
        continue

    # CHECK INITIAL CONDITIONS

    c1 = gap <= gap_max
    c2 = len(sta_list) >= nst
    c3 = rms <= rms_max
    c4 = her <= her_max
    c5 = zer <= zer_max

    if all([c1, c2, c3, c4, c5]):
        
        w_gap = check(1 - (gap/gap_max)  / gap_max_wt)
        w_nst = check(nst/len(sta_list)  * num_sta_min_wt)
        w_rms = check(1 - (rms/rms_max)  / rms_max_wt)
        w_her = check(1 - (her/her_max)  / her_max_wt)
        w_zer = check(1 - (zer/zer_max)  / zer_max_wt)
        w_smd = check(1 - (smd/sta_mdis) / sta_mdis_wt)
        
        w_tot = w_gap + w_nst + w_rms + w_her + w_zer + w_smd
        
        w_mean+=w_tot

        for cell in sel_dic:

            xmin = sel_dic[cell]['geo']['x'][0]
            xmax = sel_dic[cell]['geo']['x'][1]

            ymin = sel_dic[cell]['geo']['y'][0]
            ymax = sel_dic[cell]['geo']['y'][1]

            zmin = sel_dic[cell]['geo']['z'][0]
            zmax = sel_dic[cell]['geo']['z'][1]

            cx = xmin <= lon <= xmax
            cy = ymin <= lat <= ymax
            cz = zmin <= dep <= zmax
            
            if all([cx, cy, cz]) and evt not in sel_dic[cell]['evt_list']:

                sel_dic[cell]['evt_list'][evt] = w_tot

                break

print ''

w_mean  = w_mean / float(len(data))
sel_evt = []
            
for cell in sel_dic:

    if sel_dic[cell]['evt_list']:

        x = sel_dic[cell]['evt_list']
        x = sorted(x.items(), key=operator.itemgetter(1), reverse=True)

        for i in range(max_nevt):

            try:

                if x[i][1] >= (w_mean * (evt_wt/100.)):

                    sel_evt.append(x[i][0])

            except IndexError:

                pass
        
for i in sel_evt:

    x_fin.append(data[i]['HEADER']['L1']['Lon'])
    y_fin.append(data[i]['HEADER']['L1']['Lat'])
    z_fin.append(-data[i]['HEADER']['L1']['Dep'])

#___________________PLOT RESULTS

init_plotting_isi(18,9)
plt.rcParams['axes.labelsize'] = 6
plt.rcParams['axes.titlesize'] = 6

#_____AX1
    
ax1 = plt.subplot(2,2,1)
[i.set_linewidth(0.6) for i in ax1.spines.itervalues()]
ax1.tick_params(labelsize=6)
ax1.locator_params(axis='x',nbins=5)
ax1.locator_params(axis='y',nbins=5)

ax1.plot(x_ini, y_ini, marker='o', color='r', ms=2, ls='',  mew=.1, mec='k')

x,y = get_cells(sel_dic,2)

for i,j in zip(x,y): ax1.plot(i, j, ls=':', dashes=(5, 10), lw=.1, color='gray')

ax1.set_xlim(lon_min - .05*lon_len, lon_max + .05*lon_len)
ax1.set_ylim(lat_min - .05*lat_len, lat_max + .05*lat_len)
ax1.set_title('Initial (All Events, %d)'%(len(x_ini)))
ax1.locator_params(nbins=6)

#_____AX2

ax2 = plt.subplot(2,2,2)
[i.set_linewidth(0.6) for i in ax2.spines.itervalues()]
ax2.tick_params(labelsize=6)
ax2.locator_params(axis='x',nbins=5)
ax2.locator_params(axis='y',nbins=5)

ax2.plot(x_fin, y_fin, marker='o', color='r',  ms=2, ls='', mew=.1, mec='k')

for i,j in zip(x,y): ax2.plot(i, j, ls=':', dashes=(5, 10), lw=.1, color='gray')

ax2.set_xlim(lon_min - .05*lon_len, lon_max + .05*lon_len)
ax2.set_ylim(lat_min - .05*lat_len, lat_max + .05*lat_len)
ax2.set_title('Final (Selected Events, %d)'%(len(x_fin)))
ax2.locator_params(nbins=6)

#_____AX3

ax3 = plt.subplot(2,2,3,projection='3d')
[i.set_linewidth(0.1) for i in ax3.spines.itervalues()]
ax3.tick_params(labelsize=6, pad=.1)

ax3.plot(x_ini, y_ini, z_ini, marker='o', color='r', ms=2, ls='', mew=.1, mec='k')
ax3.locator_params(nbins=5)

x,y,z = get_cells(sel_dic,3)

for i,j,k in zip(x,y,z): ax3.plot(i, j, k, ls=':', dashes=(5, 10), lw=.1, color='b', alpha=.5)

ax3.set_xlim(lon_min - .05*lon_len, lon_max + .05*lon_len)
ax3.set_ylim(lat_min - .05*lat_len, lat_max + .05*lat_len)
ax3.set_zlim(-par_dic['DEP_MAX'], 0.0)

ax3.xaxis._axinfo["grid"].update({"linewidth":.5, "color" : "gray", "linestyle":"--"})
ax3.yaxis._axinfo["grid"].update({"linewidth":.5, "color" : "gray", "linestyle":"--"})
ax3.zaxis._axinfo["grid"].update({"linewidth":.5, "color" : "gray", "linestyle":"--"})

ax3.w_xaxis.set_pane_color((1.0, 1.0, 0.0, 0.5))
ax3.w_yaxis.set_pane_color((1.0, 1.0, 0.0, 0.6))
ax3.w_zaxis.set_pane_color((1.0, 1.0, 0.0, 0.7))

#_____AX4

ax4 = plt.subplot(2,2,4,projection='3d')
[i.set_linewidth(0.6) for i in ax4.spines.itervalues()]
ax4.tick_params(labelsize=6, pad=1)

ax4.plot(x_fin, y_fin, z_fin, marker='o', color='r', ms=2, ls='', mew=.1, mec='k')
ax4.locator_params(nbins=5)

x,y,z = get_cells(sel_dic,3)

for i,j,k in zip(x,y,z): ax4.plot(i, j, k, ls=':', dashes=(5, 10), lw=.1, color='b', alpha=.5)

ax4.xaxis._axinfo["grid"].update({"linewidth":.5, "color" : "gray", "linestyle":"--"})
ax4.yaxis._axinfo["grid"].update({"linewidth":.5, "color" : "gray", "linestyle":"--"})
ax4.zaxis._axinfo["grid"].update({"linewidth":.5, "color" : "gray", "linestyle":"--"})

ax4.w_xaxis.set_pane_color((1.0, 1.0, 0.0, 0.5))
ax4.w_yaxis.set_pane_color((1.0, 1.0, 0.0, 0.6))
ax4.w_zaxis.set_pane_color((1.0, 1.0, 0.0, 0.7))

plt.tight_layout()
plt.savefig('select3d.tiff')

#___________________WRITE SELECTED EVENTS

for c in sel_evt:

    out_dic[c] = data[c]
            
Write_Nordic(inp=out_dic, output='select.out')

#___________________PRINT SUMMARY

print '+++ Cell size in km^3                      :', int(cell_size)
print '+++ Horizontal Cell size in km^2           :', int(par_dic['H_CELL_SIZE'])
print '+++ Number of cells in X, Y and Z          : %d, %d, %d'%(len(X)-1, len(Y)-1, len(Z)-1)
print '+++ Total number of cells                  :', num_cell
print '+++ Max Number of event per allowded cell  :', max_nevt
print '+++ Total Number of event (initial)        :', len(data)
print '+++ Total Number of checked event          :', len(data)-tot_evt_chk
print '+++ Total Number of selected event (final) :', len(x_fin)
print '+++ Output file with selected event(s)     : select.out\n'
