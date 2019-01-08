from numpy import loadtxt
from datetime import datetime as dt
from numpy import mean, sqrt
import os
import sys
sys.path.append(os.path.join('utility'))

if not os.path.exists('select2D.par'):

    print '\n+++ Warning!\n'
    print '- No "select2D.par" file was found in current directory!'
    print '- Default values will be used and you can change them in the next run.\n'

    f = open('select2D.par','w')
    f.write("""# set paramaeters for selecting events
#
#
# Note:
#      If RMV_PHA_DIS paramaeter is used,
#      one more  relocation IS NEEDED to 
#      perform  other  filters,  simply
#      because the location and time may
#      differ  from  the  original  one.
#
###################
START_TIME         = 1900-01-01T00:00:00        # Start time of search window [1996-01-01T00:00:00].
END_TIME           = 2100-01-01T00:00:00        # End time of search window [2020-01-01T00:00:00].
FPS                = None                       # Fault Plane Solution like FOCMEC,PINV, ... .
M_TYPE             = None                       # Magnitude Type(s) like ml,mb, ... .
DIS_ID             = None                       # Distance ID(s) like L,R,D.
EV_ID              = None                       # Event ID(s) like 20040627202532.
MAG_LIM            = 0.0  ,  9.99               # Magnitude Limits like [2.0, 5.5].
LAT_LIM            = 00.00, 90.00               # Latitude Limits like [30.00, 39.00].
LON_LIM            = 00.00, 90.00               # Longitude Limits like [45.00, 55.00].
DEP_LIM            = 00.00, 999.9               # Depth Limits like [0.0, 999.9].
RMS_LIM            = 00.00, 99.99               # RMS limits like [0.0, 99.9].
NOBS               = 0    , 999                 # Number of Stations Limits like [0, 99].
HOR_ERR_LIM        = 00.00, 999.9               # Horizontal Errors Limits like [0.0, 999.9].
DEP_ERR_LIM        = 00.00, 999.9               # Depth Errors Depth Limits like [0.0, 999.9].
MIN_N_POL          = 0                          # Minimum Number of Polarities like [0].
HYPO_ANC           = None                       # Hypocenter Agencies like [0].
MAG_AGC            = None                       # Magnitude Agencies like [None].
GAP_LIM            = 0, 360                     # Gap range  like [0, 360].
MIN_DIS            = 999.9                      # Minimum distance between station to event (at least one must exist!) like [1000.0].
MEAN_DIS           = 999.9                      # Mean of distances for all station-event pairs like [3000.0].
RMV_PHA_DIS        = Pn, 0, 4                   # Change the weigth of desired phase within the specefic distance like [Pn,200,4].
S_WEIGHT           = 2                          # Set minimum S-Phase weigth to this value. type None to disable this item.                       
RMV_FIXED_D        = True                       # Remove flag (F) for fixed depth [True] or not [False].""")
    
    f.close()
    select_par  = loadtxt('select2D.par',comments='#',dtype=str,delimiter='=')
   
else:

    select_par  = loadtxt('select2D.par',comments='#',dtype=str,delimiter='=')

select_out  = open('select.out','w')
select_dict = {}

for i in select_par:
    
    select_dict[i[0].strip()] = i[1].strip()

inp_file = raw_input('\nInput NORDIC file name:\n\n')



event_list     = [[]]
event_surfed   = 0
event_bad_hdrE = 0
event_passed   = 0
event_nopassed = 0
event_noloc    = 0
event_nomag    = 0
c              = 0

with open(inp_file) as f:

    for l in f:

        if not l.strip():
            
            event_list.append([])
            c+=1
            event_surfed+=1

        if l.strip():
            
            event_list[c].append(l)

event_list = filter(None,event_list)


# re-weigth bad phases

tmp = []

for i in event_list:

    tmp.append([])

    for l in i:

        if l.strip() and (l[79:80] == ' ' or l[79:80] == '4') and \
           l[10:12] == select_dict['RMV_PHA_DIS'].split(',')[0].strip() and \
           l[70:75].strip() and float(l[70:75].strip()) <= float(select_dict['RMV_PHA_DIS'].split(',')[1].strip()):

            l = l[:14]+l[14:15].replace(l[14:15], select_dict['RMV_PHA_DIS'].split(',')[2].strip())+l[15:]

        if l.strip() and (l[79:80] == ' ' or l[79:80] == '4') and l[10:11] == 'S' and select_dict['S_WEIGHT'] != 'None' and float('0'+l[14:15]) <= float(select_dict['S_WEIGHT']):

            l = l[:14]+'%d'%(float(select_dict['S_WEIGHT']))+l[15:]


        if select_dict['RMV_FIXED_D']=='True' and l.strip() and l[79:80]=='1' and l[43:44]=='F':

            l = l[:43]+' '+l[44:]

        tmp[-1].append(l)

event_list = tmp

for i in event_list:

    flag_linetype1 = True

    y      = i[0][1:5].strip()
    m      = '%02d'%int(i[0][6:8].strip())
    d      = '%02d'%int(i[0][8:10].strip())
    H      = '%02d'%int(i[0][11:13].strip())
    M      = '%02d'%int(i[0][13:15].strip())
    S      = '%02d'%int(i[0][16:18].strip())
    F      = i[0][18:20].strip()
    
    if float(S) >= 60.0:
        
        M = '%02d'%(int(M)+1)
        S = '%02d'%(int(float(S)-60.0))

    if float(M) >= 60:

        M = '%02d'%(int(float(M) - 60))
        H = '%02d'%(int(float(H) + 1))

    if float(H) >= 24:

        H = '%02d'%(int(float(H) - 24))
        d = '%02d'%(int(float(d) + 1))
        
    T      = '-'.join([y,m,d])+'T'+':'.join([H,M,S])+F
    ORIGT  = dt.strptime(T,'%Y-%m-%dT%H:%M:%S.%f')
    DIS_ID = i[0][21:22].strip()
    LAT    = float(i[0][24:30])
    LON    = float(i[0][31:38])
    DEP    = float(i[0][38:43])
    LAGC   = i[0][45:48].strip()
    NOBS   = int(i[0][49:51].strip())
    RMS    = float(i[0][51:55])
    MAG    = float(i[0][55:59] if i[0][55:59].strip() else 0.0)
    DIST   = [float(l[70:75]) for l in i if l.strip() and l[10].upper()=="P" and l[70:75].strip() and l[79:80] == ' ' or l[79:80] == '4']
    try:
        MIN_D  = min(DIST)
    except: continue
    MN_D   = mean(DIST)
    FOC    = [l[70:76].strip() for l in i if l[79:80] == 'F']
    GAP    = [float(l[5:8]) for l in i if l[79:80] == 'E']
    OTER   = [float(l[14:20]) for l in i if l[79:80] == 'E']
    LATER  = [float(l[24:30]) for l in i if l[79:80] == 'E']
    LONER  = [float(l[32:38]) for l in i if l[79:80] == 'E']
    DEPER  = [float(l[38:43]) for l in i if l[79:80] == 'E']
    POL    = [l[16:17].strip() for l in i if l.strip() and l[79:80] == ' ' or l[79:80] == '4']

    STRT  = dt.strptime(select_dict['START_TIME'],'%Y-%m-%dT%H:%M:%S')
    ENDT  = dt.strptime(select_dict['END_TIME'],'%Y-%m-%dT%H:%M:%S')
    if select_dict['MAG_LIM']=='None':
        MINM = 0.0
        MAXM = 9.9
    else:
        MINM  = float(select_dict['MAG_LIM'].split(',')[0])
        MAXM  = float(select_dict['MAG_LIM'].split(',')[1])
    MILA  = float(select_dict['LAT_LIM'].split(',')[0])
    MALA  = float(select_dict['LAT_LIM'].split(',')[1])
    MILO  = float(select_dict['LON_LIM'].split(',')[0])
    MALO  = float(select_dict['LON_LIM'].split(',')[1])
    MIDP  = float(select_dict['DEP_LIM'].split(',')[0])
    MADP  = float(select_dict['DEP_LIM'].split(',')[1])
    MIRM  = float(select_dict['RMS_LIM'].split(',')[0])
    MARM  = float(select_dict['RMS_LIM'].split(',')[1])
    MIOB  = float(select_dict['NOBS'].split(',')[0])
    MAOB  = float(select_dict['NOBS'].split(',')[1])
    MIHE  = float(select_dict['HOR_ERR_LIM'].split(',')[0])
    MAHE  = float(select_dict['HOR_ERR_LIM'].split(',')[1])
    MIDE  = float(select_dict['DEP_ERR_LIM'].split(',')[0])
    MADE  = float(select_dict['DEP_ERR_LIM'].split(',')[1])
    MIPO  = float(select_dict['MIN_N_POL'])
    MIGAP = float(select_dict['GAP_LIM'].split(',')[0])
    MAGAP = float(select_dict['GAP_LIM'].split(',')[1])
    MIND  = float(select_dict['MIN_DIS'])
    MND   = float(select_dict['MEAN_DIS'])

    c        = 0
    Pn_list  = []
    num_pol  = 0
    dis_list = []
    flag_LT1 = True
    flag_LTE = True
    
    for l in i:

        if l[79:80] == '1' and flag_LT1:

            flag_LT1 = False

            try:

                chk_OT = bool(STRT <= ORIGT <= ENDT)

                if not l[55:59].strip():

                    MAG = 0.0
                    event_nomag+=1

                else:

                    MAG = float(l[55:59])
                    
                chk_MA = bool(MINM <= MAG <= MAXM)
                chk_LA = bool(MILA <= float(l[23:30]) <= MALA)
                chk_LO = bool(MILO <= float(l[30:38]) <= MALO)
                chk_DP = bool(MIDP <= float(l[38:43]) <= MADP)
                chk_RM = bool(MIRM <= float(l[51:55]) <= MARM)
                chk_OB = bool(MIOB <= float(l[48:51]) <= MAOB)

            except ValueError:

                chk_OT = False
                chk_MA = False
                chk_LA = False
                chk_LO = False
                chk_DP = False
                chk_RM = False
                chk_OB = False
                event_bad_hdr1+=1
                flag_linetype1 = False              

        if l[79:80] == 'E' and flag_linetype1 and flag_LTE:

            flag_LTE = False

            try:

                chk_GAP = bool(MIGAP <= float(l[5:8])   <= MAGAP)
                chk_HE  = bool(MIHE  <= sqrt(float(l[32:38])**2 + float(l[24:30])**2) <= MAHE)
                chk_DE  = bool(MIDE  <= float(l[38:43]) <= MADE)

            except ValueError:

                chk_GAP = False
                chk_HE  = False
                chk_DE  = False
                event_bad_hdrE+=1
                

        if l.strip() and l[16:17].strip() and (l[79:80] == ' ' or l[79:80] == '4'):

            chk_pol = l[16:17].strip()
            
            if chk_pol:

                num_pol+=1

        c+=1


    chk_MIN_DIS  =  bool(MIN_D <= MIND)
    chk_MEAN_DIS =  bool(MN_D <= MND)


    chk_POL  = bool(MIPO <= num_pol)
    chk_list = [chk_OT, chk_MA, chk_LA,
                chk_LO, chk_HE, chk_DP,
                chk_RM, chk_OB, chk_GAP,
                chk_DE, chk_MIN_DIS,
                chk_MEAN_DIS, chk_POL]

    if all(chk_list):

        event_passed+=1

        for l in i:

            select_out.write(l)

        select_out.write('\n')

    else:

        event_nopassed+=1


    
select_out.close()

#___________________ Print log
print ''
print 'Total number of event(s)                                :  ',event_surfed
print 'Total number of event(s) with bad header in line type E :  ',event_bad_hdrE
print 'Total number of event(s) passed all conditions          :  ',event_passed
print 'Total number of event(s) NOT passed all conditions      :  ',event_nopassed
print 'Total number of event(s) with no magnitude reported!    :  ',event_nomag
