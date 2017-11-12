import os
import sys
from shutil import rmtree, move
from numpy import random, genfromtxt, pi
from datetime import timedelta as td
from datetime import datetime as dt
from random import sample
sys.path.append(os.path.join('utility'))
from PyNordicRW import Read_Nordic
from PyNordic2CNV import Nordic2CNV

from random import normalvariate
from numpy import pi , arange, genfromtxt
from shutil import move
import os

#_________ READ INPUT PARAMETERS


if not os.path.exists('addnoise.par'):

    print '\n\n+++ No "addnoise.par" file was found.\n'
    print '+++ Creating new parameter file...\n'
    print '+++ Done.'

    anf = open('addnoise.par','w')

    anf.write('#\n')
    anf.write('# Noise Generation parameter file.\n')
    anf.write('# Set the following parameters to make noise/error\n')
    anf.write('# to travel times/epicenters.\n')
    anf.write('#\n')
    anf.write('##################################\n')
    anf.write('#\n')
    anf.write('TIME_ERR_MODE                     = FIX          # Time error mode. Type "FIX" for fixed or "GAU" for gassian error.\n')
    anf.write('TIME_ERR_P                        = 0.0          # Time error magnitude for "P" phases.\n')
    anf.write('TIME_ERR_S                        = 0.0          # Time error magnitude for "S" phases.\n')
    anf.write('FIXED_ERR_M                       = C            # Error addition mode when TIME_ERR_MODE=FIX. type "R" for random and "C" for constante.\n')
    anf.write('TIME_ERR_TO_STA                   = ALL          # Apply Time error to "ALL" stations or using station list file (e.g. selsta.dat [each station per line]).\n')
    anf.write('LOC_ERR_MODE                      = FIX          # Location error mode. Type "FIX" for fixed or "GAU" for gassian error.\n')
    anf.write('LOC_ERR_LAT                       = 5.0          # Location error magnitude for latitude.\n')
    anf.write('LOC_ERR_LON                       = 5.0          # Location error magnitude for longitude.\n')
    anf.write('LOC_ERR_DEP                       = 5.0          # Location error magnitude for depth.\n')
    anf.write('ADD_LOC_ERR_TO                    = ALL          # Add Location error to LAT,LON,DEP [ALL] or just to one parameter per event [OPE].\n')
    anf.write('LOC_ERR_TO_EVT                    = ALL          # Apply location error to "ALL" events or choose "N" events randomly.\n')
    anf.write('#\n')
    anf.close()

par     = genfromtxt('addnoise.par', comments='#', dtype=None, delimiter='=')
par_dic = {}

for _ in par:

    try:

        par_dic[_[0].strip()] = float(_[1].strip())
        
    except ValueError:
        
        par_dic[_[0].strip()] = _[1].strip()

ans = raw_input('\n+++ Select input file format: type "C" for "CNV" or "N" for "NORDIC".\n\n')

if ans.upper() == 'C':

    fa = raw_input('\n+++ Enter project name [must exist in "figs" directory]:\n\n')
    fa = os.path.join('..','figs',fa,'fin_hyp.cnv')

if ans.upper() == 'N':

    fa = raw_input('\n+++Enter NORDIC file name:\n\n')

time_err_mod    = par_dic['TIME_ERR_MODE']
time_err_p      = par_dic['TIME_ERR_P']
time_err_s      = par_dic['TIME_ERR_S']
fixed_err_m     = par_dic['FIXED_ERR_M']
time_err_to_sta = par_dic['TIME_ERR_TO_STA']
loc_err_mod     = par_dic['LOC_ERR_MODE']
loc_err_lat     = par_dic['LOC_ERR_LAT']
loc_err_lon     = par_dic['LOC_ERR_LON']
loc_err_dep     = par_dic['LOC_ERR_DEP']
add_loc_err_to  = par_dic['ADD_LOC_ERR_TO']
loc_err_to_evt  = par_dic['LOC_ERR_TO_EVT']


#_________ FUNCTIONS


def k2d(kilometer, radius=6371):

    return kilometer / (2.0 * radius * pi / 360.0)


def d2k(degrees, radius=6371):

    return degrees * (2.0 * radius * pi / 360.0)

def read_selsta(inp='selsta.dat'):

    if time_err_to_sta == 'ALL':

        return 'ALL'

    else:

        if os.path.exists(inp):

            sta_list = genfromtxt(inp, dtype=str)

            return sta_list

        else:
            
            print '\n+++ No %s file was found!'%(inp)
            sys.exit(0)
                
if ans.upper() == 'N':

    nor_out = 'select_noisy.out'


    def add_time_err(line, error_mode, error_mag, fixed_err_mode):

        time = line[18:28]
        H    = int(time[0:2])
        M    = int(time[2:4])
        S    = int(float(time[4:]))
        MS   = int(1e6*(float(time[4:]) - int(float(time[4:]))))

        if S >= 60:

            S = S - 60
            M+=1

        if M >= 60:

            M = M - 60
            H+=1

        if H >= 24:

            H = H - 24
            print '+++ Found Hour > 24 ! => Skipping phase.'
            
        if error_mode == 'GAU':

            error = random.normal(0.0,error_mag)

        elif error_mode == 'FIX':

            if fixed_err_mode == 'C':

                error = error_mag

            if fixed_err_mode == 'R':
                
                rnd   = random.normal()
                error = error_mag * rnd/abs(rnd)

        time = dt(2000,1,1,H,M,S,MS)
        time = time + td(0,error)
        time = '%02d%02d %5.2f'%(time.hour,time.minute,time.second+time.microsecond*1e-6)

        return time



    def add_loc_err(loc, error_mode, error_mag):

        if error_mode == 'GAU':

            error = random.normal(0.0, error_mag)

        elif error_mode == 'FIX':

            rnd   = random.normal()
            error = error_mag * rnd/abs(rnd)

        loc+=error

        if loc <= 0.0:

            loc+=2.0*abs(error)
        
        return loc
    

    #_________ ADD TIME ERROR


    sta_list      = read_selsta(inp=time_err_to_sta)
    output_nordic = open('tmp.out','w')


    with open(fa) as f:

        for l in f:

            if type(sta_list) == type('ALL'):
                
                if l.split() and (l[79:80] == ' ' or l[79:80] == '4') and l[10:11] == 'P':

                    time = l[18:28]
                    time = add_time_err(l, time_err_mod, time_err_p, fixed_err_m)
                    l    = l[:18]+time+l[28:]

                if l.split() and (l[79:80] == ' ' or l[79:80] == '4') and l[10:11] == 'S':

                    time = l[18:28]
                    time = add_time_err(l, time_err_mod, time_err_s, fixed_err_m)
                    l    = l[:18]+time+l[28:]
                    
                output_nordic.write(l)


            else:

                if l.split() and (l[79:80] == ' ' or l[79:80] == '4') and l[1:5].strip() in sta_list and l[10:11] == 'P':

                    time = l[18:28]
                    time = add_time_err(l, time_err_mod, time_err_p, fixed_err_m)
                    l    = l[:18]+time+l[28:]

                if l.split() and (l[79:80] == ' ' or l[79:80] == '4') and l[1:5].strip() in sta_list and l[10:11] == 'S':

                    time = l[18:28]
                    time = add_time_err(l, time_err_mod, time_err_s, fixed_err_m)
                    l    = l[:18]+time+l[28:]

                output_nordic.write(l)                

    output_nordic.close()



    #_________ ADD LOCATION ERROR

    neq           = 0
    output_nordic = open('select_noisy.out','w')


    with open('tmp.out') as f:

        for l in f:

            if l[79:80] == '1' and l[23:30].strip() and l[30:38].strip() and l[38:43].strip():

                neq+=1

    if type(loc_err_to_evt) != str:
        
        loc_err_to_evt = int(loc_err_to_evt)

        sel_eq = sample(range(neq), loc_err_to_evt)
       

    with open('tmp.out') as f:

        c = 0

        for l in f:

            if loc_err_to_evt == 'ALL':

                if l[79:80] == '1' and l[23:30].strip() and l[30:38].strip() and l[38:43].strip():

                    lat = float(l[23:30])
                    lon = float(l[30:38])
                    dep = float(l[38:43])

                    if add_loc_err_to == 'ALL':

                        lat = add_loc_err(lat, loc_err_mod, k2d(loc_err_lat))
                        lon = add_loc_err(lon, loc_err_mod, k2d(loc_err_lon))
                        dep = add_loc_err(dep, loc_err_mod, loc_err_dep)

                    if add_loc_err_to == 'OPE':

                        p = sample(['LAT', 'LON', 'DEP'], 1)[0]

                        if p == 'LAT':

                            lat = add_loc_err(lat, loc_err_mod, k2d(loc_err_lat))

                        elif p == 'LON':

                            lon = add_loc_err(lon, loc_err_mod, k2d(loc_err_lon))

                        elif p == 'DEP':

                            dep = add_loc_err(dep, loc_err_mod, loc_err_dep)

                    l   = l[:23]+'%7.3f'%(lat)+'%8.3f'%(lon)+'%5.1f'%(dep)+l[43:]

                output_nordic.write(l)


            else:

                if c in sel_eq and l[79:80] == '1' and l[23:30].strip() and l[30:38].strip() and l[38:43].strip():

                    lat = float(l[23:30])
                    lon = float(l[30:38])
                    dep = float(l[38:43])

                    lat = add_loc_err(lat, loc_err_mod, k2d(loc_err_lat))
                    lon = add_loc_err(lon, loc_err_mod, k2d(loc_err_lon))
                    dep = add_loc_err(dep, loc_err_mod, loc_err_dep)

                    l   = l[:23]+'%7.3f'%(lat)+'%8.3f'%(lon)+'%5.1f'%(dep)+l[43:]

                output_nordic.write(l)

            if l[79:80] == '1' and l[23:30].strip() and l[30:38].strip() and l[38:43].strip():

                c+=1

    output_nordic.close()

    #_________ CONVERT TO CNV AND COPY TO "velinp" DIRECTORY 

    print ''

    d  = Read_Nordic('select_noisy.out')
    ce = Nordic2CNV(d)

    move('data.cnv',os.path.join('..','velinp','noisy.cnv'))
    os.remove('tmp.out')



if ans.upper() == 'C':

    Lat_err = loc_err_lat
    Lon_err = loc_err_lon
    Dep_err = loc_err_dep
    P_err   = time_err_p
    S_err   = time_err_s
    mode    = time_err_mod
    sta_flg = time_err_to_sta
    
    nor_out = 'No'

    #_________ FUNCTIONS


    def add_time_err(t, v, mode):

        if mode == 'FIX':

            t+=v

        if mode == 'GAU':

            t+=normalvariate(t,v)

        return t

        
    def add_loc_err(a):

        c   = normalvariate(0,1)
        c   = abs(c)/c
        
        return a * c

    output = open('noisy.cnv','w')

    if os.path.exists(sta_flg):

        sta_lst = genfromtxt('selsta.dat', dtype=None)
        
    with open(fa) as f:

        for l in f:

            if (l[25:26] == 'N' or l[25:26] == 'S') and  (l[35:36] == 'E' or l[35:36] == 'W'):

                lat = float(l[17:25])
                lon = float(l[26:35])
                dep = float(l[36:43])

                if add_loc_err_to == 'ALL':

                    lat+=add_loc_err(k2d(Lat_err))
                    lon+=add_loc_err(k2d(Lon_err))
                    dep+=add_loc_err(Dep_err)


                if add_loc_err_to == 'OPE':

                    p = sample(['LAT', 'LON', 'DEP'], 1)[0]

                    if p == 'LAT':

                        lat+=add_loc_err(k2d(Lat_err))

                    elif p == 'LON':

                        lon+=add_loc_err(k2d(Lon_err))

                    elif p == 'DEP':

                        dep+=add_loc_err(Dep_err)

                if dep < 0.0:

                    dep+=2*Dep_err
                    

                l = l[:17]+'%8.4f'%(lat)+l[25:26]+'%9.4f'%(lon)+l[35:36]+'%7.2f'%(dep)+l[43:]

            elif l.strip():

                for c in arange(4,72,12, dtype=int):


                    if sta_flg == 'ALL':

                        try:

                            if l[c:c+1] == 'P':

                                t = '%6.2f'%(add_time_err(float(l[c+2:c+8]), P_err, 'FIX'))
                                l = l[:c+2]+t+l[c+8:]

                            if l[c:c+1] == 'S':

                                t = '%6.2f'%(add_time_err(float(l[c+2:c+8]), S_err, 'FIX'))
                                l = l[:c+2]+t+l[c+8:]


                        except ValueError:

                            pass

                    else:
                        
                        try:

                            if l[c-4:c] in sta_lst and l[c:c+1] == 'P':

                                t = '%6.2f'%(add_time_err(float(l[c+2:c+8]), P_err, 'FIX'))
                                l = l[:c+2]+t+l[c+8:]

                            if l[c-4:c] in sta_lst and l[c:c+1] == 'S':

                                t = '%6.2f'%(add_time_err(float(l[c+2:c+8]), S_err, 'FIX'))
                                l = l[:c+2]+t+l[c+8:]


                        except ValueError:

                            pass
                        
            output.write(l)


    output.close()
                
    move('noisy.cnv',os.path.join('..','velinp','noisy.cnv'))

#_________ PRINT SUMMARY REPORT



print ''
print '+++ Input file to add noise/errors                  : ', fa
print '+++ Output NORDIC file included noise               : ', nor_out
print '+++ Time Error mode                                 : ', time_err_mod
print '+++ Time Error added to P phase arrivals (sec)      : ', time_err_p
print '+++ Time Error added to S phase arrivals (sec)      : ', time_err_s
print '+++ Fixed Error mode (when TIME_ERR_MODE=FIX)       : ', fixed_err_m
print '+++ Station file used or select ALL stations        : ', time_err_to_sta
print '+++ Location Error mode                             : ', loc_err_mod
print '+++ Latitude Error added to event epicenter (km)    : ', loc_err_lat
print '+++ Longitude Error added to event epicenter (km)   : ', loc_err_lon
print '+++ Depth Error added to event Depth (km)           : ', loc_err_dep
print '+++ Method used for adding location error           : ', add_loc_err_to
print '+++ Number of selected events to add location error : ', loc_err_to_evt

print '\n+++ File "noisy.cnv" was copied to "velinp" directory.'
