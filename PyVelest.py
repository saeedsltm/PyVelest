import os, glob
from numpy import loadtxt, diff, ones, random, array, hstack, mean, pi, where, median, degrees, ones_like
from numpy import arange, sqrt, linalg, arccos, arctan2, rad2deg, zeros_like, abs, linspace, savez, load
from scipy import zeros, shape,  concatenate
import pylab as plt
from matplotlib import ticker, gridspec
from shutil import copy, rmtree, move
from LatLon import lat_lon as ll
from osgeo import gdal
from os import system, remove, path
import sys
sys.path.append(path.join('tools','utility'))
from datetime import datetime as dt
from PyNordicRW import Read_Nordic
from matplotlib.colors import LightSource
from matplotlib.cbook import get_sample_data
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Ellipse
import tarfile
from collections import OrderedDict
import operator
from mpl_toolkits.axes_grid1 import make_axes_locatable


plt.rc('font',family='Times New Roman')
plt.style.use('seaborn-paper')

#___________________ set global settings

def init_plotting():

    plt.rcParams['figure.figsize']       = (26, 14)
    plt.rcParams['figure.dpi']           = 75
    plt.rcParams['font.size']            = 16
    plt.rcParams['font.family']          = 'Times New Roman'
    plt.rcParams['axes.labelsize']       = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize']       = 1.5*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize']      = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize']      = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize']      = plt.rcParams['font.size']
    plt.rcParams['xtick.major.size']     = 3
    plt.rcParams['xtick.minor.size']     = 3
    plt.rcParams['xtick.major.width']    = 1
    plt.rcParams['xtick.minor.width']    = 1
    plt.rcParams['ytick.major.size']     = 3
    plt.rcParams['ytick.minor.size']     = 3
    plt.rcParams['ytick.major.width']    = 1
    plt.rcParams['ytick.minor.width']    = 1
    plt.rcParams['legend.frameon']       = True
    plt.rcParams['legend.shadow']        = True
    plt.rcParams['legend.loc']           = 'lower left'
    plt.rcParams['legend.numpoints']     = 1
    plt.rcParams['legend.scatterpoints'] = 1
    plt.rcParams['axes.linewidth']       = 1
    plt.rcParams['savefig.dpi']          = 200
    plt.rcParams['xtick.minor.visible']  = 'False'
    plt.rcParams['ytick.minor.visible']  = 'False'
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')

init_plotting()


#___________________ some usseful functions

def k2d(kilometer, radius=6371):

    return kilometer / (2.0 * radius * pi / 360.0)

def d2k(degrees, radius=6371):

    return degrees * (2.0 * radius * pi / 360.0)

def reject_outliers(data, m = 2.):
    
    d             = abs(data - median(data))
    mdev          = median(d)
    s             = d/mdev if mdev else 0.
    inx           = where(s<m)
    good_data     = data[inx]
    good_data_ind = where(s<m)
    inx           = where(s>m)
    bad_data      = data[inx]
    bad_data_ind  = where(s>m)
    
    return good_data,good_data_ind[0],bad_data,bad_data_ind[0]
     
#___________________ Main class

    
class main():

    """
    Description:
    A class for running velest in automatic mode (synthetic
    model generation, control file prepration, analysing and
    plot results).

    """

    def __init__(self):

        self.project_nm = raw_input("\n+++ Enter Project Name: <Please use the name of velocity model you're using>\n\n")

        self.input_dic = self.read_par()

    #*-*-*-*-*- Cheeking if topography data exists
    #
    # Extract tar.gz file may take times

    def extract_topo(self):

        if self.topo_flag == 'True':

            here = os.getcwd()
            os.chdir(os.path.join('gdb','db','srtm'))

            if not os.path.exists('IRN.tif'):

                print "+++ You're running 'PyVelest' for the first time on your system!"
                print "+++ You'll never see this message next time."
                print "\n--> Extracting topography data , Please wait ...\n"
                
                tar = tarfile.open('IRN.tif.tar.gz', "r:gz")
                tar.extractall()
                tar.close()

            os.chdir(here)
    
    #*-*-*-*-*- Remove old files and do empy all folders
    #
    # Remove files from velout.
    # Remove tmp folder if exist.

    def clean(self):

        if os.path.exists(os.path.join('tmp',self.project_nm)):

            rmtree(os.path.join('tmp',self.project_nm))

        if os.path.exists(os.path.join('velout',self.project_nm)):

            rmtree(os.path.join('velout',self.project_nm))


    def write_cmn(self, input_dic, output='velest.cmn', project_nm='test'):

        output     = open(output,'w')
        now        = dt.now()

        header1 = "CONTROL-FILE FOR PROGRAM  V E L E S T  (%s)"%(now.strftime('%Y-%m-%d'))
        header1 = header1.center(70,'*')
        header2 = project_nm

        # write Headers
        
        output.write(header1+'\n')
        output.write(header2+'\n')

        # write Main lines

        d   = input_dic

        for _ in d:

            if d[_] == "None":

                d[_] = '    '

        l1  = '***  olat       olon   icoordsystem      zshift   itrial ztrial    ised'
        l2  = '*** neqs   nshot   rotate'
        l3  = '*** isingle   iresolcalc'
        l4  = '*** dmax    itopo    zmin     veladj    zadj   lowveloclay'
        l5  = '*** nsp    swtfac   vpvs       nmod'
        l6  = '***   othet   xythet    zthet    vthet   stathet'
        l7  = '*** nsinv   nshcor   nshfix     iuseelev    iusestacorr'
        l8  = '*** iturbo    icnvout   istaout   ismpout'
        l9  = '*** irayout   idrvout   ialeout   idspout   irflout   irfrout   iresout'
        l10 = '*** delmin   ittmax   invertratio'
        l11 = '*** Modelfile'
        l12 = '*** Stationfile'
        l13 = '*** Seismofile'
        l14 = '*** File with region names'
        l15 = '*** File with region coordinates'
        l16 = '*** File #1 with topo data'
        l17 = '*** File #2 with topo data'
        l18 = '*** File with Earthquake data'
        l19 = '*** File with Shot data'
        l20 = '*** Main print output file'
        l21 = '*** File with single event locations'
        l22 = '*** File with final hypocenters in *.cnv format'
        l23 = '*** File with new station corrections'
        l24 = '*** File with summary cards (e.g. for plotting)'
        l25 = '*** File with raypoints'
        l26 = '*** File with derivatives'
        l27 = '*** File with ALEs'
        l28 = '*** File with Dirichlet spreads'
        l29 = '*** File with reflection points'
        l30 = '*** File with refraction points'
        l31 = '*** File with residuals'
        l32 = '******* END OF THE CONTROL-FILE FOR PROGRAM  V E L E S T  *******'

        for l in [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10]:

            output.write(l+'\n')
            l = l.split()
            l = l[1:]
            l = "    ".join([d[v] for v in l ])
            output.write(l+'\n')

        output.write(l11+'\n')
        output.write(os.path.join('velinp',self.velmod_name)+'\n')
        output.write(l12+'\n')
        output.write(d['Stationfile']+'\n')
        output.write(l13+'\n')
        output.write(d['Seismofile']+'\n')
        output.write(l14+'\n')
        output.write(d['RegionNamesfile']+'\n')
        output.write(l15+'\n')
        output.write(d['Region_coorfile']+'\n')
        output.write(l16+'\n')
        output.write(d['Topo_data_1file']+'\n')
        output.write(l17+'\n')
        output.write(d['Topo_data_2file']+'\n')
        output.write(l18+'\n')
        output.write(d['EQ_datafile']+'\n')
        output.write(l19+'\n')
        output.write(d['Shot_datafile']+'\n')
        output.write(l20+'\n')
        output.write(d['Main_printoutfile']+'\n')
        output.write(l21+'\n')
        output.write(d['Single_ev_locfile']+'\n')
        output.write(l22+'\n')
        output.write(d['Final_hypfile']+'\n')
        output.write(l23+'\n')
        output.write(d['Station_corrfile']+'\n')
        output.write(l24+'\n')
        output.write(d['Summary_cardsfile']+'\n')
        output.write(l25+'\n')
        output.write(d['Raypointsfile']+'\n')
        output.write(l26+'\n')
        output.write(d['Derivativesfile']+'\n')
        output.write(l27+'\n')
        output.write(d['ALEsfile']+'\n')
        output.write(l28+'\n')
        output.write(d['Dirichlet_sprfile']+'\n')
        output.write(l29+'\n')
        output.write(d['Reflection_pntfile']+'\n')
        output.write(l30+'\n')
        output.write(d['Refraction_pntfile']+'\n')
        output.write(l31+'\n')
        output.write(d['Residualsfile']+'\n')
        output.write(l32+'\n')

        output.close()

        return project_nm



    #*-*-*-*-*- read VELEST input parameters into a dictionary
    #
    # Given the parameter file [./par/par.dat], read all VELEST
    # required parameters.

    def read_par(self,par_file=os.path.join('par','par.dat')):   

        print "\n+++ Read VELEST input parameters ..."

        self.input_par = loadtxt(par_file,comments='#',dtype=str,delimiter=':')
        self.input_dic = {}

        for i in self.input_par:

            self.input_dic[i[0].strip()] = i[1].strip()


        # Read topo.dat
        
        parfile = loadtxt(path.join('par','topo.dat'),skiprows=3,comments='#',dtype=str,delimiter='=')
        par_dic = {}

        for _ in parfile:
            
            par_dic[_[0].strip()] = _[1].strip()

        self.lon_min, self.lon_max = float(par_dic['LON_MIN']), float(par_dic['LON_MAX'])
        self.lat_min, self.lat_max = float(par_dic['LAT_MIN']), float(par_dic['LAT_MAX'])
        self.dep_min, self.dep_max = float(par_dic['DEP_MIN']), float(par_dic['DEP_MAX']) 
        self.mag_min, self.mag_max = float(par_dic['MAG_MIN']), float(par_dic['MAG_MAX'])
        self.mag_flag              = par_dic['MAG_FLG']
        self.max_her               = float(par_dic['MAX_HER'])
        self.bin_her               = float(par_dic['BIN_HER'])
        self.max_der               = float(par_dic['MAX_DER'])
        self.bin_der               = float(par_dic['BIN_DER'])
        self.max_rms               = float(par_dic['MAX_RMS'])
        self.max_tdf               = float(par_dic['TIMEDIF'])
        self.max_ldf               = float(par_dic['LOCDIF'])
        self.gmtcpt                = par_dic['GMT_CPT']
        self.pltcpt                = eval('plt.cm.'+par_dic['PLT_CPT'])
        self.cptflag               = par_dic['GMTPLTF']
        self.srtm_tif              = par_dic['SRTMTIF']
        self.res_x                 = par_dic['RES_P_X']
        self.res_y                 = par_dic['RES_P_Y']
        self.topo_flag             = par_dic['TOPOFLG']
        self.nstd                  = float(par_dic['NSTD'])
        self.sc_len                = int(float(par_dic['SC_LEN']))
        self.bat_dep               = par_dic['BAT_DEP']
        self.plt_flt               = par_dic['PLT_FLT']
        self.ot_adj                = float(par_dic['OT_ADJ'])
        self.hr_adj                = float(par_dic['HR_ADJ'])
        self.dp_adj                = float(par_dic['DP_ADJ'])

        # Read velocity models in velinp directory and let user to choose wich one must be used

        inp_syntvel      = os.path.join('par','syntvel.dat')
        inp              = loadtxt(inp_syntvel,dtype=str,comments='#',delimiter=':')
        self.velmod_name = inp[16][1].strip()

        return self.input_dic
    


    #*-*-*-*-*- write new VELEST parameter to velest.cmn file
    #
    # Write VELEST parameters into VELEST's cnv
    # file [./velinp/velest.cmn]
    
    def write_par(self):

        self.input_dic = self.read_par()

        self.project_nm = self.write_cmn(input_dic=self.input_dic, output='velest.cmn',project_nm=self.project_nm)

        print "+++ Write VELEST input parameters ..."


    #*-*-*-*-*- read cnv data file into a dictionary
    #
    # Read a cnv file format and parse it into a dictionary
    # New dictionary include, origin time, lat,lon,depth,mag
    # azimuthal gap and rms for each event followded by a nested
    # dictionary includes all phase arrivals.

    def read_cnv(self,inpcnv=os.path.join('velinp','data.cnv')):

        print "+++ Read phase data ..."

        cnv_dic = {}

        with open(inpcnv) as f:

            for l in f:

                if l.strip():

                    chk_hdr = [h.isdigit() or h==' ' for h in l[0:5]]

                    if all(chk_hdr):

                        yer  = int('%004d'%float(l[0:2]))
                        mon  = int('%01d'%float(l[2:4]))
                        day  = int('%1d'%float(l[4:6]))
                        H    = int('%1d'%float(l[7:9]))
                        M    = int('%1d'%float(l[9:11]))
                        sec  = float(l[12:17])
                        msec = int((sec - int(sec))*1e6)
                        sec  = int(sec)

                        if sec >= 60:

                            sec = sec - 60
                            M+=1

                        if M >= 60:

                            M = M - 60
                            H+=1

                        if H >= 24:

                            H = H - 24
                            day+=1

                        if yer == 0:

                            yer = 1900
                                
                        try:

                            ot   = dt(yer,mon,day,H,M,sec,msec)

                        except ValueError:

                            print '+++ ==> Error encounterd on the following line:\n',l
                            
                        lat  = float(l[18:26].strip()[:-1])
                        lon  = float(l[27:36].strip()[:-1])
                        dep  = float(l[37:43].strip())
                        mag  = float(l[45:50].strip())
                        gap  = float(l[54:57].strip())
                        rms  = float(l[62:67].strip())
                        eqid = ot

                        cnv_dic[eqid] = {'lat':lat,
                                         'lon':lon,
                                         'dep':dep,
                                         'mag':mag,
                                         'gap':gap,
                                         'rms':rms,
                                         'arr':{'sta':[],
                                                'pha':[],
                                                'wt':[],
                                                'tim':[]}
                                         }

                    elif not all(chk_hdr):

                        cc = 0

                        for _ in range(6):

                            try:

                                cnv_dic[eqid]['arr']['sta'].append(l[0+cc:4+cc].strip())
                                cnv_dic[eqid]['arr']['pha'].append(l[4+cc:5+cc].strip())
                                cnv_dic[eqid]['arr']['wt'].append(int(l[5+cc:6+cc].strip()))
                                cnv_dic[eqid]['arr']['tim'].append(float(l[7+cc:12+cc]))

                            except ValueError:

                                pass
                            
                            cc+=12 # length of each phase info

                        cnv_dic[eqid]['arr']['sta'] = filter(None,cnv_dic[eqid]['arr']['sta'])
                        cnv_dic[eqid]['arr']['pha'] = filter(None,cnv_dic[eqid]['arr']['pha'])
                        cnv_dic[eqid]['arr']['tim'] = filter(None,cnv_dic[eqid]['arr']['tim'])

        return cnv_dic
                            


    #*-*-*-*-*- read velocity model and station files
    #
    # Read VELEST velocity file [./velinp/model.mod] and
    # station file [./velinp/station.sta] and parse it into
    # a dictionary.

    def read_VelSta(self,inp_vel=os.path.join('velinp','model.mod'), inp_sta=os.path.join('velinp','station.sta')):

        inp_vel=os.path.join('velinp',self.velmod_name)

        print "+++ Read velocity and station information ..."

        sta_dic = {}

        with open(inp_sta) as f:

            for l in f:

                if '(a4,f7.4,a1,1x,f8.4,a1,1x,i4,1x,i1,1x,i3,1x,f5.2,2x,f5.2)' not in l and l.strip():

                    nm  = l[0:4]
                    lat = float(l[4:11])
                    lon = float(l[13:21])
                    elv = float(l[22:27])
                    rfn = int(l[30:33].strip())

                    sta_dic[nm] = {'lat':lat, 'lon':lon, 'elv':elv, 'rfn':rfn}
                    

        with open(inp_vel) as f:

            vel_dic = {}

            p_vel,p_dep,p_vdamp = [], [], []
            s_vel,s_dep,s_vdamp = [], [], []

            start_p_flag = False
            start_s_flag = False

            for l in f:

                if "P-VELOCITY MODEL" in l:

                    start_p_flag = True
                
                if "S-VELOCITY MODEL" in l:

                    start_s_flag = True
                    start_p_flag = False

                if start_p_flag and len(l.split()) >= 3:

                    p_vel.append(float(l.split()[0]))
                    p_dep.append(float(l.split()[1]))
                    p_vdamp.append(float(l.split()[2]))

                    vel_dic['P'] = {'p_vel':p_vel, 'p_dep':p_dep, 'p_dmp':p_vdamp}

                if start_s_flag and len(l.split()) >= 3:

                    s_vel.append(float(l.split()[0]))
                    s_dep.append(float(l.split()[1]))
                    s_vdamp.append(float(l.split()[2]))

                    vel_dic['S'] = {'s_vel':s_vel, 's_dep':s_dep, 's_dmp':s_vdamp}

        return vel_dic,sta_dic        



    #*-*-*-*-*- Make synthetic velocity models
    #
    # Read parameter file's ./par/syntvel.dat and make
    # synthetic velocity models

    def mk_syn_model(self,inp_syntvel=os.path.join('par','syntvel.dat')):

        print "+++ Make synthetic velocity models ..."

        inp       = loadtxt(inp_syntvel,dtype=str,comments='#',delimiter=':')
        
        vel_dic,_ = self.read_VelSta(inp_vel=os.path.join('velinp',self.velmod_name),inp_sta=self.input_dic['Stationfile'])
        pvel      = vel_dic['P']['p_vel']
        pdep      = vel_dic['P']['p_dep']
        pdmp      = vel_dic['P']['p_dmp']

        try:

            sdmp      = vel_dic['S']['s_dmp']

        except KeyError:

           sdmp = ones(len(pdmp)) 
            
        self.vel_mode      = inp[0][1].strip()
        self.dep_mode      = inp[1][1].strip()
        self.vel_min       = float(inp[2][1])
        self.vel_max       = float(inp[3][1])
        self.dep_min       = float(inp[4][1])
        self.dep_max       = float(inp[5][1])
        self.vel_dev       = float(inp[6][1])
        self.dep_dev       = float(inp[7][1])
        self.nmod          = int(inp[8][1].strip())
        self.vpvs          = float(inp[9][1])
        self.ors           = float(inp[10][1])
        self.grd           = int(inp[11][1].strip())
        self.bad_evnt_herr = int(inp[12][1].strip())
        self.bad_evnt_zerr = int(inp[13][1].strip())
        self.bad_evnt_rms  = float(inp[14][1].strip())
        self.res_max_rw    = float(inp[15][1].strip())
        self.velmod_name    = inp[16][1].strip()
       
        
        if not os.path.exists(os.path.join('tmp',self.project_nm,'synthvel')):

            os.makedirs(os.path.join('tmp',self.project_nm,'synthvel'))

        else:

            rmtree(os.path.join('tmp',self.project_nm,'synthvel'))
            os.makedirs(os.path.join('tmp',self.project_nm,'synthvel'))
            
        
        output      = os.path.join('tmp',self.project_nm,'synthvel')
        orig_velmod = array([vel_dic['P']['p_vel'],vel_dic['P']['p_dep'],vel_dic['P']['p_dmp']]).T

        # generate random model with no gradient
        
        if self.grd == 0:
            
            if self.vel_mode == "u":

                vel = array([random.uniform(i-self.vel_dev*(1.0/j),i+self.vel_dev*(1.0/j),self.nmod) for i,j in zip(pvel,pdmp)]).T

            elif self.vel_mode == "n":

                vel = array([random.normal(i,self.vel_dev*(1.0/j),self.nmod) for i,j in zip(pvel,pdmp)]).T


            if self.dep_mode == "u":

                dep = array([random.uniform(i-self.dep_dev*(1.0/j),i+self.dep_dev*(1.0/j),self.nmod) for i,j in zip(pdep,pdmp)]).T

            elif self.dep_mode == "n":

                dep = array([random.normal(i,self.dep_dev*(1.0/j),self.nmod) for i,j in zip(pdep,pdmp)]).T


        # generate random model with positive gradient

        elif self.grd == 1:


            vel = []
            dep = []
            
            if self.vel_mode == "u":

                for _ in range(self.nmod):

                    vel.append(sorted(random.uniform(i-self.vel_dev*(1.0/j),i+self.vel_dev*(1.0/j)) for i,j in zip(orig_velmod[:,0],pdmp)))

            elif self.vel_mode == "n":

                for _ in range(self.nmod):

                    vel.append(sorted(random.normal(i,self.vel_dev*(1.0/j))  for i,j in zip(orig_velmod[:,0],pdmp)))


            if self.dep_mode == "u":

                for _ in range(self.nmod):

                    dep.append(sorted(random.uniform(i-self.dep_dev*(1.0/j),i+self.dep_dev*(1.0/j))  for i,j in zip(orig_velmod[:,1],pdmp)))

            elif self.dep_mode == "n":

                for _ in range(self.nmod):

                    dep.append(sorted(random.normal(i,self.dep_dev*(1.0/j))  for i,j in zip(orig_velmod[:,1],pdmp)))

            vel = array(vel)
            dep = array(dep)            


        dep_org     = [i[1] for i in orig_velmod]          
        vel_org     = [i[0] for i in orig_velmod]
        dep_org[0]  = self.dep_min

        for d in range(len(dep)):
            
            dep[d][0] = self.dep_min      

            
        nl = len(dep[0]) # number of layers for each random model
        

        # Write to a new model file
        
        
        for i in range(self.nmod):
            
            fid = file(os.path.join(output,"model_"+str(i+1)+".mod"),'w')
            fid.write('%s%d%s\n'%("Model: model_",i+1,".mod"))
            fid.write('%s%d %s\n'%(" ",len(dep[0]),"         vel,depth,vdamp,phase(f5.2,5x,f7.2,2x,f7.3,3x,a1)"))
            
            for j in range(len(dep[0])):
                
                if j==0:
                    
                    fid.write('%5.2f %11.2f %8.3f %26s\n'%(vel[i][j],dep[i][j],pdmp[j],"P-VELOCITY MODEL"))
                    
                else:
                    
                    fid.write('%5.2f %11.2f %8.3f\n'%(vel[i][j],dep[i][j],1.0))
                    
            fid.write('%s%s\n'%(" ",len(dep[0])))
            
            for j in range(len(dep[0])):
                
                if j==0:
                    
                    fid.write('%5.2f %11.2f %8.3f %26s\n'%(vel[i][j]/self.vpvs,dep[i][j],pdmp[j],"S-VELOCITY MODEL"))
                    
                else:
                    
                    fid.write('%5.2f %11.2f %8.3f\n'%(vel[i][j]/self.vpvs,dep[i][j],1.0))
                    
            fid.close()

        # Plot synthetic models

        plt.figure()

        
        for i in range(self.nmod):

            y  = [-1*d for d in dep[i]]
            yy = [-1*d for d in dep_org]
            
            x  = list(hstack([[m,n] for m,n in zip(vel[i],vel[i])]))
            xx = list(hstack([[m,n] for m,n in zip(vel_org,vel_org)]))
            y  = list(hstack([[m,n] for m,n in zip(y,y)]))
            yy = list(hstack([[m,n] for m,n in zip(yy,yy)]))
            
            y.pop(0)
            yy.pop(0)

            y.append(-1*self.dep_max)        
            yy.append(-1*self.dep_max)

            if i == 0:
                
                plt.plot(x, y, 'r-', linewidth=2, label='Random')
                
            else:
                
                plt.plot(x, y, 'r-', linewidth=2)
                plt.ylim(ymin=-self.dep_max,ymax=-self.dep_min)
                plt.title("Synthetic Models")
                plt.xlabel("Velocity [km/s]")
                plt.ylabel("Depth [km]")

        velmod_org = [xx,yy]
        
        plt.plot(xx, yy, 'k-', linewidth=3,label='Reference')
        plt.grid()
        plt.legend()
        plt.xlim(self.vel_min, self.vel_max)
        plt.tight_layout()
        plt.savefig(os.path.join('figs',self.project_nm,"synt_models_"+self.project_nm))
        plt.show()
        plt.close()


    #*-*-*-*-*- Run VELEST program
    #
    # Run VELEST for each model
    # 

    def run_velest(self):

        print "+++ Run VELEST ..."

        if not os.path.exists(os.path.join('velout',self.project_nm)):

            os.makedirs(os.path.join('velout',self.project_nm))

        if os.path.exists(os.path.join('tmp',self.project_nm,'synthvel')):

            self.synthvel_files = glob.glob(os.path.join('tmp',self.project_nm,'synthvel','*.mod'))
            self.synthvel_files = sorted(self.synthvel_files, key=lambda x: int(x.split(os.sep)[-1].split('_')[1].split('.')[0]))

            print '   --- %d synthetic model(s) found ...'%(len(self.synthvel_files))

        else:

            self.synthvel_files = []
            self.nmod           = 0

        if self.synthvel_files:

            copy('velest.cmn',os.path.join('tmp',self.project_nm))
    
            for m,c in zip(self.synthvel_files,range(len(self.synthvel_files))):

                tmp_cmn = open(os.path.join('tmp',self.project_nm,'tmp.cmn'),'w')
                
                with open(os.path.join('tmp',self.project_nm,'velest.cmn')) as f:

                    for line in f:

                        if os.path.join('velinp',self.velmod_name) in line:

                            tmp_cmn.write(line.replace(os.path.join('velinp',self.velmod_name), m))

                        elif os.path.join('velout','velest.out') in line:

                            tmp_cmn.write(line.replace(os.path.join('velout','velest.out'),
                                                       os.path.join('velout',self.project_nm,'velest'+str(c+1)+'.out')))

                        elif os.path.join('velout','final_loc.cnv') in line:

                            tmp_cmn.write(line.replace(os.path.join('velout','final_loc.cnv'),
                                                       os.path.join('velout',self.project_nm,'final_loc'+str(c+1)+'.cnv')))

                        elif os.path.join('velout','stations_corr.sta') in line:
                            tmp_cmn.write(line.replace(os.path.join('velout','stations_corr.sta'),
                                                       os.path.join('velout',self.project_nm,'stations_corr'+str(c+1)+'.sta')))

                        else:

                            tmp_cmn.write(line)

                tmp_cmn.close()
                
                move(os.path.join('tmp',self.project_nm,'tmp.cmn'),'velest.cmn')

                print '   --- Run for synthetic model:',c+1

                os.system('velest > /dev/null')

            os.remove(os.path.join('tmp',self.project_nm,'velest.cmn'))

            for _ in [os.path.join('velout','velest.out'),
                      os.path.join('velout','final_loc.cnv'),
                      os.path.join('velout','stations_corr.sta')]:

                if os.path.exists(_):

                    os.remove(_)


        else:

            print '   --- No synthetic model was found ...'
            print '   --- Run for initial velocity model ...'
            os.system('velest > /dev/null')

        

    #*-*-*-*-*- Analyse VELEST results 
    #
    # Analyse VELEST results and plot 
    # them.
    
    def analyse_velest_res(self):

        print "+++ Analyse VELEST outputs ..."

        inp_syntvel    = os.path.join('par','syntvel.dat')
        inp            = loadtxt(inp_syntvel,dtype=str,comments='#',delimiter=':')
        self.vel_mode  = inp[0][1].strip()
        self.dep_mode  = inp[1][1].strip()
        self.vel_min   = float(inp[2][1])
        self.vel_max   = float(inp[3][1])
        self.dep_min   = float(inp[4][1])
        self.dep_max   = float(inp[5][1])
        self.vel_dev   = float(inp[6][1])
        self.dep_dev   = float(inp[7][1])
        self.nmod      = int(inp[8][1].strip())
        self.vpvs      = float(inp[9][1])
        self.nit       = int(self.input_dic['ittmax'])

        # if multiple synthetic velocity model has been defined
        
        if self.synthvel_files:

            self.rms_res     = [[] for sm in self.synthvel_files]
            self.avg_abs_adj = [[] for sm in self.synthvel_files]
            self.vel_adj     = [[] for sm in self.synthvel_files]
            self.svel        = {'P':[[] for sm in self.synthvel_files],
                                'S':[[] for sm in self.synthvel_files]}
            self.sdep        = {'P':[[] for sm in self.synthvel_files],
                                'S':[[] for sm in self.synthvel_files]}
            self.stacor      = {'P':[[] for sm in self.synthvel_files],
                                'S':[[] for sm in self.synthvel_files]}

            s_vel, s_dep     = [], [] # if nsp = 1, just one model will be  produced

            # read each VELEST output to extract inversion results
            
            for sm in range(self.nmod):
                

                with open(os.path.join('velout',self.project_nm,'velest%d.out'%(sm+1))) as f:

                    self.p_vel_adj_flag = False
                    self.s_vel_adj_flag = False
                    vel_flag            = False
                    iter_counter        = 0

                    for l in f:

                        if 'DATVAR=' in l and iter_counter <= self.nit:

                            self.rms_res[sm].append(float(l[68:].strip()))
                            iter_counter+=1

                        if 'A V E R A G E   of ABSOLUTE ADJUSTMENTS' in l and iter_counter-1 <= self.nit:

                            self.avg_abs_adj[sm].append([float(i) for i in l[42:].strip().split()])

                        if  'Velocity model   1' in l and iter_counter-1 <= self.nit:

                            self.p_vel_adj_flag = True
                            tmp_vel             = []
                            p_vel, p_dep        = [], []

                        if  'Velocity model   2' in l and iter_counter-1 <= self.nit:

                            self.s_vel_adj_flag = True
                            s_vel, s_dep        = [], []
                            
                        if self.p_vel_adj_flag and 'Velocity model   1' not in l and len(l.split())==3:

                            tmp_vel.append((float(l.split()[2]),float(l.split()[1])))
                            p_vel.append(float(l.split()[0]))
                            p_dep.append(float(l.split()[2]))

                        if self.s_vel_adj_flag and 'Velocity model   2' not in l and len(l.split())==3:

                            s_vel.append(float(l.split()[0]))
                            s_dep.append(float(l.split()[2]))
                                         
                        if not l.strip() and self.p_vel_adj_flag:
                            
                            self.p_vel_adj_flag = False
                            self.vel_adj[sm].append(tmp_vel)

                        if not l.strip() and self.s_vel_adj_flag:
                            
                            self.s_vel_adj_flag = False

                self.svel['P'][sm].append(p_vel)
                self.sdep['P'][sm].append(p_dep)
                self.svel['S'][sm].append(s_vel)
                self.sdep['S'][sm].append(s_dep)    

            self.p_outlier_models = []
            self.s_outlier_models = []
            self.p_good_models    = []
            self.s_good_models    = []
            vp                    = array([_[0] for _ in self.svel['P']])
            vs                    = array([_[0] for _ in self.svel['S']])

            for _ in range(len(vp[0])):

                good_model, good_model_ind, bad_model, bad_model_ind = reject_outliers(vp[:,_],m=self.ors)
                self.p_outlier_models.extend(list(bad_model_ind))

            for _ in range(len(vs[0])):

                good_model, good_model_ind, bad_model, bad_model_ind = reject_outliers(vs[:,_],m=self.ors)
                self.s_outlier_models.extend(list(bad_model_ind))

            self.p_outlier_models = list(set(self.p_outlier_models))
            self.p_good_models    = list(set(range(self.nmod)) - set(self.p_outlier_models))
            self.svel_p           = [i[0] for i in self.svel['P']]
            self.sdep_p           = [i[0] for i in self.sdep['P']]
            self.nomc             = len(self.p_good_models)

            self.s_outlier_models = list(set(self.s_outlier_models))
            self.s_good_models    = list(set(range(self.nmod)) - set(self.s_outlier_models))

            
            # plot results, Time/Location/Velocity adjusment 

            output     = os.path.join('figs','invres_'+self.project_nm)
            iter_color = ['k','r','g','b','c','y','m','purple','pink','brown'] # iteration color, maximum 10
                               
            ax1         = plt.subplot(2,3,1)
            self.minmod = {} 
            for nm in self.p_good_models:
                x      = range(1,1+self.nit)
                ax1.plot([self.avg_abs_adj[nm][_][0] for _ in range(len(self.avg_abs_adj[nm]))], marker='o',color='grey')
                self.minmod[nm] = mean([self.avg_abs_adj[nm][_][0] for _ in range(len(self.avg_abs_adj[nm]))])
            self.ind = sorted(self.minmod.items(), key=operator.itemgetter(1))[0][0]
            self.minmod = [self.avg_abs_adj[self.ind][i][0] for i in range(len(self.avg_abs_adj[self.ind]))]
            ax1.plot(self.minmod,marker='o',color='red',label='Minimum model,%d'%(self.ind))
            ax1.set_xlabel('Iteration Number [#]')
            ax1.set_ylabel('Origin time Adjustment [sec]')
            ax1.set_ylim(0.0,self.ot_adj)
            ax1.set_xlim(0,self.nit-1)
            plt.grid()
            ax1.locator_params(axis='x',nbins=5)
            ax1.locator_params(axis='y',nbins=5)
            ax1.legend(loc=1)

            ax2 = plt.subplot(2,3,2)
            self.minmod = {}
            for nm in self.p_good_models:
                x      = range(1,1+self.nit)
                ax2.plot([self.avg_abs_adj[nm][_][1] for _ in range(len(self.avg_abs_adj[nm]))], marker='o',color='grey')
                self.minmod[nm]=mean([self.avg_abs_adj[nm][_][1] for _ in range(len(self.avg_abs_adj[nm]))])
            self.ind = sorted(self.minmod.items(), key=operator.itemgetter(1))[0][0]
            self.minmod = [self.avg_abs_adj[self.ind][i][1] for i in range(len(self.avg_abs_adj[self.ind]))]
            ax2.plot(self.minmod,marker='o',color='red',label='Minimum model,%d'%(self.ind))
            ax2.set_xlabel('Iteration Number [#]')
            ax2.set_ylabel('Longitude Adjustment [km]')
            ax2.set_ylim(0.0,self.hr_adj)
            ax2.set_xlim(0,self.nit-1)
            plt.grid()
            ax2.locator_params(axis='x',nbins=5)
            ax2.locator_params(axis='y',nbins=5)
            ax2.legend(loc=1)

            ax3 = plt.subplot(2,3,3)
            self.minmod = {}
            for nm in self.p_good_models:
                x      = range(1,1+self.nit)
                ax3.plot([self.avg_abs_adj[nm][_][2] for _ in range(len(self.avg_abs_adj[nm]))], marker='o',color='grey')
                self.minmod[nm]=mean([self.avg_abs_adj[nm][_][2] for _ in range(len(self.avg_abs_adj[nm]))])
            self.ind = sorted(self.minmod.items(), key=operator.itemgetter(1))[0][0]
            self.minmod = [self.avg_abs_adj[self.ind][i][2] for i in range(len(self.avg_abs_adj[self.ind]))]
            ax3.plot(self.minmod,marker='o',color='red',label='Minimum model,%d'%(self.ind))
            ax3.set_xlabel('Iteration Number [#]')
            ax3.set_ylabel('Latitude Adjustment [km]')
            ax3.set_ylim(0.0,self.hr_adj)
            ax3.set_xlim(0,self.nit-1)
            plt.grid()
            ax3.locator_params(axis='x',nbins=5)
            ax3.locator_params(axis='y',nbins=5)
            ax3.legend(loc=1)

            ax4 = plt.subplot(2,3,4)
            self.minmod = {}
            for nm in self.p_good_models:
                x      = range(1,1+self.nit)
                ax4.plot([self.avg_abs_adj[nm][_][3] for _ in range(len(self.avg_abs_adj[nm]))], marker='o',color='grey')
                self.minmod[nm]=mean([self.avg_abs_adj[nm][_][3] for _ in range(len(self.avg_abs_adj[nm]))])
            self.ind = sorted(self.minmod.items(), key=operator.itemgetter(1))[0][0]
            self.minmod = [self.avg_abs_adj[self.ind][i][3] for i in range(len(self.avg_abs_adj[self.ind]))]
            ax4.plot(self.minmod,marker='o',color='red',label='Minimum model,%d'%(self.ind))
            ax4.set_xlabel('Iteration Number [#]')
            ax4.set_ylabel('Depth Adjustment [km]')
            ax4.set_ylim(0.0,self.dp_adj)
            ax4.set_xlim(0,self.nit-1)
            plt.grid()
            ax4.locator_params(axis='x',nbins=5)
            ax4.locator_params(axis='y',nbins=5)
            ax4.legend(loc=1)

            ax5 = plt.subplot(2,3,5)
            self.minmod = {}
            self.test = {}
            for nm in self.p_good_models:
                x      = range(1,1+self.nit)
                ax5.plot(self.rms_res[nm],marker='o',color='grey')
                self.minmod[nm]=mean(self.rms_res[nm])
            self.ind = sorted(self.minmod.items(), key=operator.itemgetter(1))[0][0]
            self.minmod = self.rms_res[self.ind]
            ax5.plot(self.minmod,marker='o',color='red',label='Minimum model,%d'%(self.ind))
            ax5.set_xlabel('Iteration Number [#]')
            ax5.set_ylabel('RMS-Residual [sec]')
            plt.grid()
            ax5.locator_params(axis='x',nbins=5)
            ax5.locator_params(axis='y',nbins=5)
            ax5.legend(loc=1)


            ax6 = plt.subplot(2,3,6)

            x,y, ymax = [],[],[]

            for it,c,s in zip(self.vel_adj[self.ind],range(len(iter_color)),linspace(3,20,self.nit)[::-1]):
                
                color = iter_color[c]
                                
                for p in it:
                    
                    x.append(p[0])
                    y.append(p[1])
                    
                ymax.append(max([abs(_) for _ in y]))
                it = array(it)
                ax6.plot([it[:,1],zeros_like(it[:,1])],[it[:,0],it[:,0]],
                         markevery=it[:,0].size,lw=s,c=color,label='Iter: '+str(c+1))

            ax6.set_xlabel('Velocity Adjusment [km/s]')
            ax6.set_ylabel('Depth [km]')
            ax6.set_xlim(-max(ymax),max(ymax))
            ax6.locator_params(axis='x',nbins=5)
            ax6.locator_params(axis='y',nbins=5)
            ax6.set_ylim(ax6.get_ylim()[::-1])
            plt.grid()
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            lgnd=plt.legend(by_label.values(), by_label.keys(),fontsize=12)
            for legobj in lgnd.legendHandles:
                legobj.set_linewidth(5.0)
            ax6.axvline(lw=2,color='k',linestyle='--')

            plt.tight_layout()
            plt.savefig(os.path.join('figs',self.project_nm,'velest_analytics_'+self.project_nm))
            plt.close()


            # plot results, final velocity, station correction beside final location using minimum model
            
            cnv_dic         = self.read_cnv(inpcnv=os.path.join('velout',self.project_nm,'final_loc%d.cnv'%(self.ind+1)))
            vel_dic,sta_dic = self.read_VelSta(inp_vel=os.path.join('velinp',self.velmod_name),inp_sta=self.input_dic['Stationfile'])

            print "+++ Plot data statistics ..."

            ovel  = [] # origin
            odep  = [] # origin
            stnm  = sta_dic.keys()
            stlat = [sta_dic[i]['lat'] for i in sta_dic]                # station latitude
            stlon = [sta_dic[i]['lon'] for i in sta_dic]                # station longitude
            stref = [sta_dic[i]['rfn'] for i in sta_dic]                # station reference
            evlat = [cnv_dic[i]['lat'] for i in cnv_dic]                # event latitude
            evlon = [cnv_dic[i]['lon'] for i in cnv_dic]                # event longitude
            evdep = [cnv_dic[i]['dep'] for i in cnv_dic]                # event depth
            evgap = [cnv_dic[i]['gap'] for i in cnv_dic]                # event gap
            evrms = [cnv_dic[i]['rms'] for i in cnv_dic]                # event rms
            evarr = [cnv_dic[i]['arr'] for i in cnv_dic]                # event arrival time
            nump  = [cnv_dic[i]['arr']['wt'].count(0)+
                     cnv_dic[i]['arr']['wt'].count(1)+
                     cnv_dic[i]['arr']['wt'].count(2)+
                     cnv_dic[i]['arr']['wt'].count(3) for i in cnv_dic] # event arrival time number             
            refid = stref.index(max(stref))

            
            for _,__ in zip(range(len(vel_dic['P']['p_vel'])),range(len(vel_dic['P']['p_dep']))):

                ovel.append(vel_dic['P']['p_vel'][_])
                ovel.append(vel_dic['P']['p_vel'][__])
                odep.append(vel_dic['P']['p_dep'][__])
                odep.append(vel_dic['P']['p_dep'][__])
                
            odep.pop(0)
            odep.append(self.dep_max)

            tmpvel_p = [[] for nm in range(self.nmod)]
            tmpdep_p = [[] for nm in range(self.nmod)]

            for nm in range(self.nmod):
                
                for _,__ in zip(self.svel_p[nm],self.sdep_p[nm]):
                    
                    tmpvel_p[nm].append(_)
                    tmpvel_p[nm].append(_)
                    tmpdep_p[nm].append(__)
                    tmpdep_p[nm].append(__)
                
                tmpdep_p[nm].pop(0)
                tmpdep_p[nm].append(self.dep_max)

                self.svel_p[nm] = tmpvel_p[nm]
                self.sdep_p[nm] = tmpdep_p[nm]

            
            # plot-No 1, velocity model

            gs = gridspec.GridSpec(4, 4)

            fig = plt.figure()

            ax1 = fig.add_subplot(gs[0:2,0:2])
            
            for nm in self.p_good_models:
                ax1.plot(self.svel_p[nm],self.sdep_p[nm],color='red',linewidth=2)
            ax1.plot(self.svel_p[0],self.sdep_p[0],color='red',linewidth=2,label='Inverted')
            ax1.plot(ovel,odep,color='black',linewidth=2.5,label='Refrence')
            ax1.plot(self.svel_p[self.ind],self.sdep_p[self.ind],color='green',linewidth=2,label='Minimum')
            ax1.plot(mean(self.svel_p,axis=0),mean(self.sdep_p,axis=0),color='blue',linewidth=1.5,label='Mean')
            xlim = ax1.get_xlim()
            ylim = ax1.get_ylim()
            xlim = (self.vel_min,self.vel_max)
            ylim = [min(odep),max(odep)+0.1*diff(ylim)[0]]
            ax1.set_xlim(xlim)
            ax1.set_ylim(ylim)
            ax1.set_ylim(ax1.get_ylim()[::-1])
            ax1.set_xlabel('Velocity [km/s]')
            ax1.set_ylabel('Depth [km]')
            plt.grid()
            ax1.legend(numpoints=1)
            ax1.locator_params(axis='x',nbins=8)
            ax1.locator_params(axis='y',nbins=5)

            # plot results, station correction

            stnm  = []
            stlat = []
            stlon = []
            stcor = []

            with open(os.path.join('velout',self.project_nm,'stations_corr%d.sta'%(self.ind+1))) as f:

                for l in f:

                    if l.strip() and '(a4,f7.4,a1,1x,f8.4,a1,1x,i4,1x,i1,1x,i3,1x,f5.2,2x,f5.2)' not in l:

                        stnm.append(l[0:4])
                        stlat.append(float(l[4:11]))
                        stlon.append(float(l[13:21]))
                        stcor.append(float(l[34:39]))

            stnm  = array(stnm)
            stlat = array(stlat)
            stlon = array(stlon)
            stcor = array(stcor)

            # Plot station corrections
            
            ax2 = fig.add_subplot(gs[0:2,2:])

            # plot some features

            self.m = Basemap(llcrnrlon=self.lon_min, llcrnrlat=self.lat_min,
                        urcrnrlon=self.lon_max, urcrnrlat=self.lat_max,fix_aspect=False)
            self.tehran = loadtxt(path.join('gdb','db','IRN_PROV','tehran.poly'))
            self.x,self.y    = self.tehran[:,0], self.tehran[:,1]
            self.m.plot(self.x,self.y, marker=None,color='r',linewidth=2.0,zorder=1)
            ax2.text(51.3890, 35.6892, '$Tehran$', clip_on=True)
            ax2.plot(51.3700, 35.6892,marker='s', markersize=5, color='k', clip_on=True)


            #____ Faults, Borders

            self.m.readshapefile(path.join('gdb','db','IRN_PROV','IRN_adm0'),'iran')
            self.m.readshapefile(path.join('gdb','db','IRN_PROV','IRN_adm1'),'prov', color='grey')

            #____ Ticks

            parallels = linspace(self.lat_min, self.lat_max, 5)
            meridians = linspace(self.lon_min, self.lon_max, 5)

            #____ Labels = [left,right,top,bottom]

            parallels = linspace(self.lat_min,self.lat_max,5)
            self.m.drawparallels(parallels,labels=[True,False,False,False],fmt='%4.1f',linewidth=0)
            meridians = linspace(self.lon_min,self.lon_max,5)
            self.m.drawmeridians(meridians,labels=[False,False,False,True],fmt='%4.1f',linewidth=0)

            #____ Scale & North sing

            self.scale_len = self.sc_len # km
            self.sc_lon_st = min(ax2.get_xlim())+0.03*diff(ax2.get_xlim())[0]
            self.sc_lon_ed = self.sc_lon_st + k2d(self.scale_len)
            self.sc_lat_st = min(ax2.get_ylim())+0.03*diff(ax2.get_ylim())[0]
            self.sc_lat_ed = self.sc_lat_st
            self.txt_x     = mean([self.sc_lon_st, self.sc_lon_ed])
            self.txt_y     = self.sc_lat_ed + abs(diff([self.sc_lat_ed,min(ax2.get_ylim())]))

            ax2.plot([self.sc_lon_st,self.sc_lon_ed],
                    [self.sc_lat_st,self.sc_lat_ed],
                    color='k',linewidth=3.0)
            ax2.plot(self.sc_lon_st,self.sc_lat_ed,
                    color='k',linewidth=3.0, marker='>',ms=10)
            ax2.plot(self.sc_lon_ed,self.sc_lat_ed,
                    color='k',linewidth=3.0, marker='<',ms=10)
            ax2.text(self.txt_x,self.txt_y,str(self.scale_len)+' km',
                    horizontalalignment='center')

            ax2.annotate('N', xy=(.030, .07), xycoords='axes fraction', fontsize=16,  ha='center', va='bottom')
            ax2.annotate('^', xy=(.031, .09), xycoords='axes fraction', fontsize=22,  ha='center', va='bottom')

            scl = 500
            c   = evrms
            im  = ax2.scatter(evlon,evlat,c=c,marker='.',s=100,edgecolor='black',cmap=plt.cm.cool,vmin=0.0, vmax=self.max_rms)
            ax2.scatter(stlon[0],stlat[0],s=0.5*scl,marker='+',color='none',linewidth=1,label='+%.1f'%(self.max_rms))
            ax2.scatter(stlon[stcor>=0],stlat[stcor>=0],s=abs(stcor[stcor>=0]*scl),marker='+',color='blue',linewidth=1)
            ax2.scatter(stlon[0],stlat[0],s=0.5*scl,marker='o',color='none',linewidth=1,label='-%.1f'%(self.max_rms))
            ax2.scatter(stlon[stcor<0] ,stlat[stcor<0], s=abs(stcor[stcor<0]*scl),marker='o',facecolors='none' ,linewidth=1,edgecolor='red')

            for x,y,t in zip (stlon,stlat,stnm):

                if self.lon_min<=x<=self.lon_max and self.lat_min<=y<=self.lat_max:

                    ax2.text(x,y,t,fontsize=10)

            plt.grid()
            ax2.set_xlim(self.lon_min,self.lon_max)
            ax2.set_ylim(self.lat_min,self.lat_max)
            ax2.locator_params(axis='x',nbins=5)
            ax2.locator_params(axis='y',nbins=5)
            leg = ax2.legend(loc=4)
            legend = ax2.get_legend()
            legend.legendHandles[0].set_color('red')
            legend.legendHandles[1].set_edgecolor('blue')
            divider      = make_axes_locatable(ax2)
            cax          = divider.append_axes("right", size="2%", pad=0.05)
            cb           = plt.colorbar(im, ax=ax2, cax=cax, label='RMS [second]')
            tick_locator = ticker.MaxNLocator(nbins=6)
            cb.locator   = tick_locator
            cb.update_ticks()

            # plot compare figure between original and minimum model

            self.orig_cnv = self.read_cnv(inpcnv=os.path.join('velinp','data.cnv'))
            self.mmod_cnv = self.read_cnv(inpcnv=os.path.join('velout',self.project_nm,'final_loc%d.cnv'%(self.ind+1)))

            olon = []
            olat = []
            odep = []
            orms = []

            mlon = []
            mlat = []
            mdep = []
            mrms = []

            lat_diff = []
            lon_diff = []
            
            for i in sorted(self.orig_cnv.keys()):

                i_lon = self.orig_cnv[i]['lon']
                i_lat = self.orig_cnv[i]['lat']
                i_dep = self.orig_cnv[i]['dep']
                i_rms = self.orig_cnv[i]['rms']

                for f in sorted(self.mmod_cnv.keys()):

                    f_lon  = self.mmod_cnv[f]['lon']
                    f_lat  = self.mmod_cnv[f]['lat']
                    f_dep  = self.mmod_cnv[f]['dep']
                    f_rms  = self.mmod_cnv[f]['rms']
                    
                    dx     = f_lon - i_lon
                    dy     = f_lat - i_lat
                    diff_h = sqrt(dx**2 + dy**2)
                    
                    diff_ot = i-f
                    diff_ot = abs(diff_ot.total_seconds())
                    
                    if diff_ot <= self.max_tdf and diff_h <= self.max_ldf:
                        
                        olon.append(i_lon)
                        olat.append(i_lat)
                        odep.append(i_dep)
                        orms.append(i_rms)

                        mlon.append(f_lon)
                        mlat.append(f_lat)
                        mdep.append(f_dep)
                        mrms.append(f_rms)

                        break

            olon = array(olon)
            olat = array(olat)
            odep = array(odep)
            orms = array(orms)

            mlon = array(mlon)
            mlat = array(mlat)
            mdep = array(mdep)
            mrms = array(mrms)
            
            ax3 = fig.add_subplot(gs[2:,0])
            
            for (x,xx) in zip(olon,mlon):
                
                if xx-x>0:

                    ax3.plot([x,x],[d2k(xx-x),0],'r-o',markevery=olon.size)

                else:

                    ax3.plot([x,x],[d2k(xx-x),0],'b-o',markevery=olon.size)
                    
            ax3.set_xlabel('Longitude [deg]')
            ax3.set_ylabel('Dislocation [km]')
            ax3.set_ylim(-max([abs(_) for _ in ax3.get_ylim()]),max([abs(_) for _ in ax3.get_ylim()]))
            ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            xlim=ax3.get_xlim()
            ax3.plot(xlim,[0,0], linewidth=2, color = 'k')
            ax3.set_xlim(xlim)
            ax3.locator_params(axis='x',nbins=5)
            ax3.locator_params(axis='y',nbins=5)
            plt.grid()

            ax4 = fig.add_subplot(gs[2:,1])

            for (y,yy) in zip(olat,mlat):

                if yy-y>0:

                    ax4.plot([y,y],[d2k(yy-y),0],'r-o',markevery=olat.size)

                else:

                    ax4.plot([y,y],[d2k(yy-y),0],'b-o',markevery=olat.size)

            ax4.set_xlabel('Latitude [deg]')
            ax4.set_ylabel('Dislocation [km]')
            ax4.set_ylim(-max([abs(_) for _ in ax4.get_ylim()]),max([abs(_) for _ in ax4.get_ylim()]))
            ax4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            xlim=ax4.get_xlim()
            ax4.plot(xlim,[0,0], linewidth=2, color = 'k')
            ax4.set_xlim(xlim)
            ax4.locator_params(axis='x',nbins=5)
            ax4.locator_params(axis='y',nbins=5)
            plt.grid()
            
            ax5 = fig.add_subplot(gs[2:,2])

            for (x,y) in zip(odep,mdep):

                if y-x>0:

                    ax5.plot([x,x],[y-x,0],'r-o',markevery=odep.size)
                else:

                    ax5.plot([x,x],[y-x,0],'b-o',markevery=odep.size)

            ax5.set_xlabel('Depth [km]')
            ax5.set_ylabel('Dislocation [km]')
            ax5.set_ylim(-max([abs(_) for _ in ax5.get_ylim()]),max([abs(_) for _ in ax5.get_ylim()]))
            ax5.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            xlim=(0,self.dep_max)
            ax5.plot(xlim,[0,0], linewidth=2, color = 'k')
            ax5.set_xlim(xlim)
            ax5.locator_params(axis='x',nbins=5)
            ax5.locator_params(axis='y',nbins=5)
            plt.grid()
            
            ax6 = fig.add_subplot(gs[2:,3])
            bns = arange(0, self.max_rms + self.max_rms/20.0, self.max_rms/20.0)
            ax6.hist(orms,color='r',label='Initial',alpha=.65,bins=bns)
            ax6.hist(mrms,color='b',label='Relocated',alpha=.65,bins=bns)
            ax6.set_xlabel('RMS [sec]')
            ax6.set_ylabel('Event [#]')
            ax6.set_xlim(0,self.max_rms)
            ax6.locator_params(axis='x',nbins=5)
            ax6.locator_params(axis='y',nbins=5)
            plt.grid()
            ax6.legend(loc=1)

            
            print '+++ Minimum model id:',self.ind+1
            plt.tight_layout()
            plt.savefig(os.path.join('figs',self.project_nm,'velest_final_'+self.project_nm))
            plt.close()


                        
        # if no synthetic velocity model has been defined (single mode)

        else:

            pass


    #*-*-*-*-*- Plot data figures
    #
    # Plot some statistical figures
    # about data.
    
    def plot_data(self):

        inp_syntvel    = os.path.join('par','syntvel.dat')
        inp            = loadtxt(inp_syntvel,dtype=str,comments='#',delimiter=':')
        self.vel_mode  = inp[0][1].strip()
        self.dep_mode  = inp[1][1].strip()
        self.vel_min   = float(inp[2][1])
        self.vel_max   = float(inp[3][1])
        self.dep_min   = float(inp[4][1])
        self.dep_max   = float(inp[5][1])
        self.vel_dev   = float(inp[6][1])
        self.dep_dev   = float(inp[7][1])
        self.nmod      = int(inp[8][1].strip())
        self.vpvs      = float(inp[9][1])  

        if not os.path.exists('figs'):

            os.mkdir('figs')

        if os.path.exists(os.path.join('figs',self.project_nm)):

            rmtree(os.path.join('figs',self.project_nm))
            os.mkdir(os.path.join('figs',self.project_nm))
    
        else:

            os.mkdir(os.path.join('figs',self.project_nm))

        cnv_dic         = self.read_cnv(self.input_dic['EQ_datafile'])
        vel_dic,sta_dic = self.read_VelSta(inp_vel=os.path.join('velinp',self.velmod_name),inp_sta=self.input_dic['Stationfile'])


        self.sta_list = {}

        with open(os.path.join('tools','STATION0.HYP')) as f:

            for l in f:

                if (l[13:14] == 'N' or l[13:14] == 'S')  and (l[22:23] == 'E' or l[22:23] == 'W'):

                    n   = l[2:6].strip()
                    lat = ll.Latitude(degree=float(l[6:8]), minute=float(l[8:13])).decimal_degree
                    lon = ll.Longitude(degree=float(l[14:17]), minute=float(l[17:22])).decimal_degree
                    
                    self.sta_list[n] = {'lon':lon,'lat':lat}
                    
        print "+++ Plot data statistics ..."

        vel   = []
        dep   = []
        stnm  = sta_dic.keys()
        stlat = [sta_dic[i]['lat'] for i in sta_dic]                # station latitude
        stlon = [sta_dic[i]['lon'] for i in sta_dic]                # station longitude
        stref = [sta_dic[i]['rfn'] for i in sta_dic]                # station reference
        evlat = [cnv_dic[i]['lat'] for i in cnv_dic]                # event latitude
        evlon = [cnv_dic[i]['lon'] for i in cnv_dic]                # event longitude
        evdep = [cnv_dic[i]['dep'] for i in cnv_dic]                # event depth
        evgap = [cnv_dic[i]['gap'] for i in cnv_dic]                # event gap
        evrms = [cnv_dic[i]['rms'] for i in cnv_dic]                # event rms
        evarr = [cnv_dic[i]['arr'] for i in cnv_dic]                # event arrival time
        nump  = [cnv_dic[i]['arr']['wt'].count(0)+
                 cnv_dic[i]['arr']['wt'].count(1)+
                 cnv_dic[i]['arr']['wt'].count(2)+
                 cnv_dic[i]['arr']['wt'].count(3) for i in cnv_dic] # event arrival time number
        refid = stref.index(max(stref))

        for _,__ in zip(range(len(vel_dic['P']['p_vel'])),range(len(vel_dic['P']['p_dep']))):

            vel.append(vel_dic['P']['p_vel'][_])
            vel.append(vel_dic['P']['p_vel'][__])
            dep.append(vel_dic['P']['p_dep'][__])
            dep.append(vel_dic['P']['p_dep'][__])
            
        dep.pop(0)
        dep.append(self.dep_max)


        # plot-No 1, velocity model

        ax1 = plt.subplot(2,2,1)
        ax1.plot(vel,dep,color='black',linewidth=2,label='Refrence')
        xlim = ax1.get_xlim()
        ylim = ax1.get_ylim()
        xlim = (self.vel_min,self.vel_max)
        ylim = [min(dep),max(dep)+0.1*diff(ylim)[0]]
        ax1.set_xlim(xlim)
        ax1.set_ylim(ylim)
        ax1.set_ylim(ax1.get_ylim()[::-1])
        ax1.set_xlabel('Velocity [km/s]')
        ax1.set_ylabel('Depth [km]')
        plt.grid()
        ax1.legend(numpoints=1)
        ax1.locator_params(axis='x',nbins=8)
        ax1.locator_params(axis='y',nbins=5)

        # plot-No 2, station-event distribution

        ax2 = plt.subplot(2,2,2)

        # plot some features

        self.m        = Basemap(llcrnrlon=self.lon_min, llcrnrlat=self.lat_min,urcrnrlon=self.lon_max, urcrnrlat=self.lat_max,
                                fix_aspect=False)
        self.tehran   = loadtxt(path.join('gdb','db','IRN_PROV','tehran.poly'))
        self.x,self.y = self.tehran[:,0], self.tehran[:,1]

        self.m.plot(self.x,self.y, marker=None,color='r',linewidth=2.0,zorder=1)
        
        ax2.text(51.3890, 35.6892, '$Tehran$')
        ax2.plot(51.3700, 35.6892,marker='s', markersize=5, color='k')

        #____ Faults, Borders

        self.m.readshapefile(path.join('gdb','db','IRN_PROV','IRN_adm0'),'iran')
        self.m.readshapefile(path.join('gdb','db','IRN_PROV','IRN_adm1'),'prov', color='grey')

        #____ Ticks on parallels and meridians

        parallels = linspace(self.lat_min, self.lat_max, 5)
        meridians = linspace(self.lon_min, self.lon_max, 5)

        #____ Labels = [left,right,top,bottom]

        parallels = linspace(self.lat_min,self.lat_max,5)
        meridians = linspace(self.lon_min,self.lon_max,5)

        self.m.drawparallels(parallels,labels=[True,False,False,False],fmt='%4.1f',linewidth=0)
        self.m.drawmeridians(meridians,labels=[False,False,False,True],fmt='%4.1f',linewidth=0)

        #____ Scale & North sing

        self.scale_len = self.sc_len # km
        self.sc_lon_st = min(ax2.get_xlim())+0.03*diff(ax2.get_xlim())[0]
        self.sc_lon_ed = self.sc_lon_st + k2d(self.scale_len)
        self.sc_lat_st = min(ax2.get_ylim())+0.03*diff(ax2.get_ylim())[0]
        self.sc_lat_ed = self.sc_lat_st
        self.txt_x     = mean([self.sc_lon_st, self.sc_lon_ed])
        self.txt_y     = self.sc_lat_ed + abs(diff([self.sc_lat_ed,min(ax2.get_ylim())]))

        ax2.plot([self.sc_lon_st,self.sc_lon_ed],
                [self.sc_lat_st,self.sc_lat_ed],
                color='k',linewidth=3.0)
        ax2.plot(self.sc_lon_st,self.sc_lat_ed,
                color='k',linewidth=3.0, marker='>',ms=10)
        ax2.plot(self.sc_lon_ed,self.sc_lat_ed,
                color='k',linewidth=3.0, marker='<',ms=10)
        ax2.text(self.txt_x,self.txt_y,str(self.scale_len)+' km',
                horizontalalignment='center')

        ax2.annotate('N', xy=(.030, .07), xycoords='axes fraction', fontsize=16,  ha='center', va='bottom')
        ax2.annotate('^', xy=(.031, .09), xycoords='axes fraction', fontsize=22,  ha='center', va='bottom')

        c   = evrms
        im  = ax2.scatter(evlon,evlat,c=c,marker='.',s=75,linewidth=None,edgecolor='black',cmap=plt.cm.cool, vmin=0.0, vmax=self.max_rms)
        ax2.scatter(stlon,stlat, s=100, marker='^', edgecolor='red', linewidth=1.5, facecolor='y',alpha=0.5)
        ax2.scatter(stlon[refid],stlat[refid], s=150, marker='^', edgecolor='g', linewidth=1.5, facecolor='r',label='Reference')

        for x,y,t in zip (stlon,stlat,stnm):
            
            if self.lon_min<=x<=self.lon_max and self.lat_min<=y<=self.lat_max:
                
                ax2.text(x,y-.03,t,fontsize=10, ha='center',bbox=dict(facecolor='w', edgecolor='k', boxstyle='round,pad=.2'))
                
        ax2.grid(True)
        ax2.locator_params(axis='x',nbins=5)
        ax2.locator_params(axis='y',nbins=5)
        ax2.legend(loc=1, fontsize=12)
        divider      = make_axes_locatable(ax2)
        cax          = divider.append_axes("right", size="2%", pad=0.05)
        cb           = plt.colorbar(im, ax=ax2, cax=cax, label='RMS [second]')
        tick_locator = ticker.MaxNLocator(nbins=6)
        cb.locator   = tick_locator
        cb.update_ticks()
        ax2.set_xlim(self.lon_min,self.lon_max)
        ax2.set_ylim(self.lat_min,self.lat_max)
        
        # plot-No 3, depth-number_of_station distribution

        ax3 = plt.subplot(2,2,3)
        c   = evrms
        im  = ax3.scatter(evdep,nump,c=c,marker='.',s=100,edgecolor='black',cmap=plt.cm.cool, vmin=0.0, vmax=self.max_rms)
        plt.grid()
        ax3.set_xlabel('Depth [km]')
        ax3.set_ylabel('Number of Phases [P,S]')
        ax3.set_xlim(0,max(evdep))
        ax3.locator_params(axis='x',nbins=5)
        ax3.locator_params(axis='y',nbins=5)
        divider      = make_axes_locatable(ax3)
        cax          = divider.append_axes("right", size="2%", pad=0.05)
        cb           = plt.colorbar(im, ax=ax3, cax=cax, label='RMS [second]')
        tick_locator = ticker.MaxNLocator(nbins=6)
        cb.locator   = tick_locator
        cb.update_ticks()

        # plot-No 4, gap-number_of_station distribution

        ax4 = plt.subplot(2,2,4)
        c   = evrms
        im  = ax4.scatter(evgap,nump,c=c,marker='.',s=100,edgecolor='black',cmap=plt.cm.cool, vmin=0.0, vmax=self.max_rms)
        plt.grid()
        ax4.set_xlabel('Azimuthal Gap [deg]')
        ax4.set_ylabel('Number of Phases [P,S]')
        ax4.set_xlim(min(evgap),max(evgap))
        ax4.locator_params(axis='x',nbins=5)
        ax4.locator_params(axis='y',nbins=5)
        divider      = make_axes_locatable(ax4)
        cax          = divider.append_axes("right", size="2%", pad=0.05)
        cb           = plt.colorbar(im, ax=ax4, cax=cax, label='RMS [second]')
        tick_locator = ticker.MaxNLocator(nbins=6)
        cb.locator   = tick_locator
        cb.update_ticks()
        
        plt.tight_layout()
        plt.savefig(os.path.join('figs',self.project_nm,'data_'+self.project_nm))
        plt.show()
        plt.close()


    ############################################## 
    #*-*-*-*-*- Make final reports
    #
    # Make final reports
    # them.
    
    def mk_report(self):

        print "+++ Preparing reports ..."

        self.repfi    = open(os.path.join('figs',self.project_nm,'report.dat'),'w')
        self.modout   = open(os.path.join('figs',self.project_nm,'model.mod'),'w')

        velocity_res  = self.svel_p
        dep_res       = self.sdep_p
        minmod_id     = self.ind
        self.min_mod  = [self.svel_p[self.ind],self.sdep_p[self.ind]]
        self.mean_mod = [mean(self.svel_p,axis=0),mean(self.sdep_p,axis=0)]

        self.repfi.write('Sumarry Report on final results\n\n')
        self.repfi.write('# 1- Minimum velocity model\n\n')
        self.repfi.write('ID number                       : %d\n'%(int(self.ind)+1))
        self.repfi.write('Number of model converged       : %d\n'%(self.nomc))
        self.repfi.write('Number of Iteration             : %d\n'%(self.nit))

        self.repfi.write('\nIteration  RMS [sec]\n')

        for i,r in zip(range(self.nit+1),self.rms_res[self.ind]):
            
            self.repfi.write('  %2d       %5.3f\n'%(i,r))
            
                         
        self.repfi.write('\nVel[km/s]   Dep[km]\n')

        self.min_mod[1][0] = 0.0 # set first depth to zero
        
        for v,d in zip(self.min_mod[0][::2],self.min_mod[1][::2]):
            
            self.repfi.write(' %4.1f      %4.1f\n'%(v,d))

        self.repfi.write('\n# 2- Mean velocity model:\n\n')

        self.repfi.write('Vel[km/s]   Dep[km]\n')

        self.mean_mod[1][0] = 0.0 # set first depth to zero
        
        for v,d in zip(self.mean_mod[0][::2],self.mean_mod[1][::2]):
            
            self.repfi.write(' %4.1f      %4.1f\n'%(v,d))
            
        self.repfi.close()

        # copy final loc and station_corr files
        
        copy(os.path.join('velout',self.project_nm,'final_loc'+str(self.ind+1)+'.cnv'),os.path.join('figs',self.project_nm,'fin_hyp.cnv'))
        copy(os.path.join('velout',self.project_nm,'stations_corr'+str(self.ind+1)+'.sta'),os.path.join('figs',self.project_nm,'sta_cor.out'))


        # write min and mean velocity model for next velest run

        self.min_mod[1][0]  = self.dep_min # set first depth to highest elevation
        self.mean_mod[1][0] = self.dep_min # set first depth to highest elevation

        # min model
        
        self.modout.write(' Model: min.mod\n')
        self.modout.write(' %2d        vel,depth,vdamp,phase(f5.2,5x,f7.2,2x,f7.3,3x,a1)\n'%(len(self.min_mod[0][::2])))

        # P-Velocity
        
        flag = True
        
        for v,d in zip(self.min_mod[0][::2],self.min_mod[1][::2]):

            if flag:

                self.modout.write(' %4.2f       %5.2f    1.000           P-VELOCITY MODEL\n'%(v,d))
                
                flag = False

                continue
                
            self.modout.write(' %4.2f       %5.2f    1.000\n'%(v,d))
            
        self.modout.write(' %2d\n'%(len(self.min_mod[0][::2])))

        # S-Velocity
        
        flag = True
        
        for v,d in zip(self.min_mod[0][::2],self.min_mod[1][::2]):
            
            if flag:

                self.modout.write(' %4.2f       %5.2f    1.000           S-VELOCITY MODEL\n'%(v/self.vpvs,d))
                
                flag = False

                continue
                
            self.modout.write(' %4.2f       %5.2f    1.000\n'%(v/self.vpvs,d))

        # now append mean model

        self.modout.write('\n')
        
        self.modout.write(' Model: mean.mod\n')
        self.modout.write(' %2d        vel,depth,vdamp,phase(f5.2,5x,f7.2,2x,f7.3,3x,a1)\n'%(len(self.mean_mod[0][::2])))

        # P-Velocity
        
        flag = True
        
        for v,d in zip(self.mean_mod[0][::2],self.mean_mod[1][::2]):

            if flag:

                self.modout.write(' %4.2f       %5.2f    1.000           P-VELOCITY MODEL\n'%(v,d))
                
                flag = False

                continue
                
            self.modout.write(' %4.2f       %5.2f    1.000\n'%(v,d))
            
        self.modout.write(' %2d\n'%(len(self.mean_mod[0][::2])))

        # S-Velocity
       
        flag = True
       
        for v,d in zip(self.mean_mod[0][::2],self.mean_mod[1][::2]):
            
            if flag:

                self.modout.write(' %4.2f       %5.2f    1.000           S-VELOCITY MODEL\n'%(v/self.vpvs,d))

                flag = False

                continue
                
            self.modout.write(' %4.2f       %5.2f    1.000\n'%(v/self.vpvs,d))

        self.modout.write('\n')
        
        self.modout.close()
        
        # locate using new velociy model and make a comparison
        
        here = os.getcwd()

        hrz_err_i = []
        dep_err_i = []
        rms_i     = []

        hrz_err_m = []
        dep_err_m = []
        rms_m     = []

        hrz_err_f = []
        dep_err_f = []
        rms_f     = []

        hrz_err_rw = []
        dep_err_rw = []
        rms_rw     = []
        
        if os.path.exists(os.path.join('tmp',self.project_nm,'report')):
            
            rmtree(os.path.join('tmp',self.project_nm,'report'))
            os.mkdir(os.path.join('tmp',self.project_nm,'report'))
            
        else:
            
            os.mkdir(os.path.join('tmp',self.project_nm,'report'))


        os.mkdir(os.path.join('tmp',self.project_nm,'report','initial'))
        os.mkdir(os.path.join('tmp',self.project_nm,'report','final'))

        # copy "select.out" and "STATION0.HYP" from "tools" directory
        
        copy(os.path.join('tools','select.out'),os.path.join('tmp',self.project_nm,'report','initial'))
        copy(os.path.join('tools','STATION0.HYP'),os.path.join('tmp',self.project_nm,'report','initial'))

        # make initial result

        hyp_par = open(os.path.join('tmp',self.project_nm,'report','initial','hyp.par'),'w')
        hyp_par.write('select.out\n')
        hyp_par.write('n\n\n')
        hyp_par.close()

        os.chdir(os.path.join('tmp',self.project_nm,'report','initial'))
        os.system('hyp < hyp.par > /dev/null')
        os.remove('gmap.cur.kml')
        os.remove('hypmag.out')

        with open('hypsum.out') as f:

            for l in f:

                if 'Date' not in l:

                    rms_i.append(float(l[62:67]))
                    hrz_err_i.append(float(l[67:72]))
                    dep_err_i.append(float(l[72:77]))
                    
        os.chdir(here)

        # make final result
        # copy "select.out" from "tools" directory
        
        copy(os.path.join('tools','select.out'),os.path.join('tmp',self.project_nm,'report','final'))

        vel        = []
        vel_flag   = True
        start_flag = False

        # make new station0 based on new velocity mode

        with open(os.path.join('figs',self.project_nm,'report.dat')) as f:

            for l in f:

                if 'Vel[km/s]   Dep[km]' in l:
                    
                    if vel_flag:

                        vel_flag = False
                        
                    else:

                        vel_flag = True
                        start_flag = True

                if start_flag and not l.strip():

                    start_flag = False
                                   
                if 'Vel[km/s]   Dep[km]' not in l and vel_flag and start_flag:
                    
                    vel.append(l)
                
        copy(os.path.join('figs',self.project_nm,'sta_cor.out'),os.path.join('tmp',self.project_nm,'report','final'))
                
        os.chdir(os.path.join('tmp',self.project_nm,'report','final'))

        sta_cor = []
        

        with open('sta_cor.out') as f:

            for l in f:

                if l.strip() and '(' not in l:

                    sta_n = l [0:4]
                    sta_lat = ll.Latitude(float(l[4:11]))
                    sta_lat = '%02d%05.2f%1s'%(sta_lat.degree,sta_lat.decimal_minute,sta_lat.get_hemisphere())
                    sta_lon = ll.Longitude(float(l[13:21]))
                    sta_lon = '%02d%05.2f%1s'%(sta_lon.degree,sta_lon.decimal_minute,sta_lon.get_hemisphere())
                    sta_elv = l[23:27]
                    sta_pt  = l[33:40]
                    sta_st  = l[41:46]
                    name = '  '+sta_n+sta_lat+' '+sta_lon+sta_elv
                    sta_cor.append(name)


        station0_hyp = open('STATION0.HYP','w')

        RESET = """RESET TEST(02)=500.0
RESET TEST(07)=-3.0                                                            
RESET TEST(08)=2.6                                                             
RESET TEST(09)=0.001                                                           
RESET TEST(11)=99.0                                                            
RESET TEST(13)=5.0                                                             
RESET TEST(34)=1.5                                                             
RESET TEST(35)=2.5                                                             
RESET TEST(36)=0.0                                                             
RESET TEST(41)=20000.0                                                         
RESET TEST(43)=5.0                                                             
RESET TEST(51)=3.6                                                             
RESET TEST(50)=1.0                                                             
RESET TEST(56)= 1.0                                                            
RESET TEST(58)= 99990.0                                                        
RESET TEST(40)=1.0                                                             
RESET TEST(60)=0.0                                                             
RESET TEST(71)=1.0                                                             
RESET TEST(75)=1.0                                                             
RESET TEST(76)=1.1725                                                            
RESET TEST(77)=0.0021                                                         
RESET TEST(78)=-2.237                                                           
RESET TEST(79)=1.0                                                             
RESET TEST(80)=3.0                                                             
RESET TEST(81)=1.0                                                             
RESET TEST(82)=1.0                                                             
RESET TEST(83)=1.0                                                             
RESET TEST(88)=1.0                                                             
RESET TEST(85)=0.1                                                             
RESET TEST(91)=0.1

"""

        station0_hyp.write(RESET)

        for s in sta_cor:

            station0_hyp.write(s+'\n')

        station0_hyp.write('\n')

        for v in vel:

            station0_hyp.write(v)
            
        station0_hyp.write('\n')
        station0_hyp.write('15.0   50. 300. %4.2f    \n'%(self.vpvs))
        station0_hyp.write('SLT\n')
        station0_hyp.close()

        copy('select.out','hyp.out')

        hyp_par = open('hyp.par','w')
        hyp_par.write('select.out\n')
        hyp_par.write('n\n\n')
        hyp_par.close()

        os.system('hyp < hyp.par > /dev/null')

        with open('hypsum.out') as f:

            for l in f:

                if 'Date' not in l:

                    rms_m.append(float(l[62:67]))
                    hrz_err_m.append(float(l[67:72]))
                    dep_err_m.append(float(l[72:77]))

        self.sta_cor = {}
        sta_wt  = {}

        with open('print.out') as f:

            start_flag = False

            for l in f:

                if ' station     p residual rms dev  no p     s residual rms dev  no s' in l:

                    start_flag = True
                

                if start_flag and not l.strip():

                    start_flag = False

                if start_flag and ' station     p residual rms dev  no p     s residual rms dev  no s' not in l:

                    self.sta_cor[l[1:5]] = {'p':'      ','s':'      '}

                    if l[15:23].strip():

                        self.sta_cor[l[1:5]]['p']='%+6.2f'%float(l[15:23])


                    if len(l.strip())==65 and l[44:52].strip():

                        self.sta_cor[l[1:5]]['s']='%+6.2f'%float(l[44:52])    
        

        tmp = open('tmp','w')

        with open('STATION0.HYP') as f:

            for l in f:

                for s in self.sta_cor:

                    if l[2:6] == s:

                        l = '  '+l[2:27]+self.sta_cor[s]['p']+' '+self.sta_cor[s]['s']+'\n'

                tmp.write(l)

        tmp.close()
        move('tmp','STATION0.HYP')

        bad_events = []

        os.system('hyp < hyp.par > /dev/null')
        with open('hypsum.out') as f:

            for l in f:

                if 'Date' not in l:

                    rms_f.append(float(l[62:67]))
                    hrz_err_f.append(float(l[67:72]))
                    dep_err_f.append(float(l[72:77]))

                    if float(l[67:72]) >= self.bad_evnt_herr or float(l[72:77]) >= self.bad_evnt_zerr or float(l[62:67]) >=self.bad_evnt_rms:

                        bad_events.append(l)

                else:

                    bad_events.append(l)

                   
        os.chdir(here)

        report = open(os.path.join('figs',self.project_nm,'report.dat'),'a')

        report.write('\n# 3- Location results using both initial and new velocity models:\n\n')
        report.write('                    * - Initial - * - Final (No ST_C) - * - Final (ST_C) - \n')
        report.write('Mean Hrz Error [km] :    %5.2f    ,      %5.2f          ,     %5.2f        :\n'%(mean(hrz_err_i),mean(hrz_err_m),mean(hrz_err_f)))
        report.write('Mean Dep Error [km] :    %5.2f    ,      %5.2f          ,     %5.2f        :\n'%(mean(dep_err_i),mean(dep_err_m),mean(dep_err_f)))
        report.write('Mean RMS      [sec] :    %5.2f    ,      %5.2f          ,     %5.2f        :\n'%(mean(rms_i),mean(rms_m),mean(rms_f)))
        report.write('\n# 4- Bad events, N = %d\n\n'%(len(bad_events)-1))
        for _ in bad_events:
            report.write(_)
        report.close()

        # relocate one more time using re-weighted phase based on the mean residual of
        # each station which has been claculated after the first relocation with station correction

        os.chdir(os.path.join('tmp',self.project_nm,'report','final'))

        sta_wt     = {}
        bad_evnt   = list(array(bad_events).copy())
        bad_events = [] 
        tmp        = open('select_rw.out','w')
                               
        # re-weight each phase on each station, base on the mean residual of that station

        for s in self.sta_cor:

            sta_wt[s] = {'p':None,'s':None}

            p_res = self.sta_cor[s]['p']
            s_res = self.sta_cor[s]['s']

            if p_res:

                p_res = self.res_max_rw
                w     = linspace(0,abs(p_res),5)
                w_dic = {(w[x],w[x+1]):x for x in range(4)}
                w_dic[(w[-1],99.9)] = 4
                sta_wt[s]['p'] = w_dic

            if s_res:

                s_res = self.res_max_rw
                w     = linspace(0,abs(s_res),5)
                w_dic = {(w[x],w[x+1]):x for x in range(4)}
                w_dic[(w[-1],99.9)] = 4
                sta_wt[s]['s'] = w_dic         
        
        with open('hyp.out') as f:

            for l in f:

                if l.strip() and (l[79:80] == ' ' or l[79:80] == '4') and l[10:11].upper() =='P' and l[63:69].strip():

                    sta_nm = l[1:5]
                    st_res = abs(float(l[63:69]))
                    try:
                        for w in sta_wt[sta_nm]['p']:
                            if w[0]<=st_res<w[1]:
                                sta_w = str(sta_wt[sta_nm]['p'][w])
                                l = l[0:14]+sta_w+l[15:]
                    except KeyError:
                        pass

                if l.strip() and (l[79:80] == ' ' or l[79:80] == '4') and l[10:11].upper() =='S' and l[63:69].strip():

                    sta_nm = l[1:5]
                    st_res = abs(float(l[63:69]))
                    try:
                        for w in sta_wt[sta_nm]['s']:
                            if w[0]<=st_res<w[1]:
                                sta_w = str(sta_wt[sta_nm]['s'][w])
                                l = l[0:14]+sta_w+l[15:]
                    except KeyError:
                        pass
                    
                tmp.write(l)

        tmp.close()

        hyp_par = open('hyp.par','w')
        hyp_par.write('select_rw.out\n')
        hyp_par.write('n\n\n')
        hyp_par.close()

        os.system('hyp < hyp.par > /dev/null')

        with open('hypsum.out') as f:

            for l in f:

                if 'Date' not in l:

                    rms_rw.append(float(l[62:67]))
                    hrz_err_rw.append(float(l[67:72]))
                    dep_err_rw.append(float(l[72:77]))

                    if float(l[67:72]) >= self.bad_evnt_herr or float(l[72:77]) >= self.bad_evnt_zerr or float(l[62:67]) >=self.bad_evnt_rms:

                        bad_events.append(l)

                else:

                    bad_events.append(l)

        os.remove('gmap.cur.kml')
        os.remove('hypmag.out')
        os.remove('hyp.par')

        os.chdir(here)

        report = open(os.path.join('figs',self.project_nm,'report.dat'),'a')

        report.write('\n# 5- Relocation results using both initial and new velocity models and rewighted phases:\n\n')
        report.write('                    * - Initial - * - Final (No ST_C) - * - Final (ST_C) - * - Final_RW (ST_C) - \n')
        report.write('Mean Hrz Error [km] :    %5.2f    ,      %5.2f          ,     %5.2f        ,     %5.2f        :\n'%(mean(hrz_err_i),
                                                                                                                          mean(hrz_err_m),
                                                                                                                          mean(hrz_err_f),
                                                                                                                          mean(hrz_err_rw)))
        report.write('Mean Dep Error [km] :    %5.2f    ,      %5.2f          ,     %5.2f        ,     %5.2f        :\n'%(mean(dep_err_i),
                                                                                                                          mean(dep_err_m),
                                                                                                                          mean(dep_err_f),
                                                                                                                          mean(dep_err_rw)))
        report.write('Mean RMS      [sec] :    %5.2f    ,      %5.2f          ,     %5.2f        ,     %5.2f        :\n'%(mean(rms_i),
                                                                                                                          mean(rms_m),
                                                                                                                          mean(rms_f),
                                                                                                                          mean(rms_rw)))
        report.write('\n# 6- Bad events, N = %d\n\n'%(len(bad_events)-1))
        for _ in bad_events:
            report.write(_)
        report.close()

        # make a new select.out with removed bad events
        
        rmd_bad_evt = open(os.path.join('figs',self.project_nm,'select.out'),'w')
        bad_evt     = open(os.path.join('figs',self.project_nm,'bad_select.out'),'w')
        event_list  = [[]]
        c           = 0

        with open(os.path.join('tools','select.out')) as f:

            flag_correct = True

            for l in f:

                if not l.strip():
                    
                    event_list.append([])
                    flag_correct = True
                    flag_LT1     = True
                    c+=1
                   
                if l.strip() and l[79:80] == '1' and not l[48:51].strip() and flag_LT1:

                    flag_LT1     = False
                    flag_correct = False
                    
                if l.strip() and flag_correct:
                    event_list[c].append(l)

        event_list = filter(None,event_list)

        for _ in event_list:

            t1 = dt.strptime(_[0][1:15],'%Y %m%d %H%M')
            flag_good = True
            flag_bad  = False

            for __ in bad_evnt:

                if ' Date' not in __:               

                    tmp = __[:11]
                    
                    if tmp[:2] == '  ':
                        
                        tmp = '00'+tmp[2:]
                        
                    else:
                        
                        tmp = '%02d'%(float(tmp[0:2]))+tmp[2:]
                        
                    t2 = dt.strptime(tmp,'%y%m%d %H%M')

                    df = t2-t1
                                 

                    if abs(df.total_seconds())<self.max_tdf:

                        flag_good = False
                        flag_bad  = True
                        
            if flag_good:

                for l in _:

                    rmd_bad_evt.write(l)
                rmd_bad_evt.write('\n')

            if flag_bad:

                for l in _:

                    bad_evt.write(l)
                bad_evt.write('\n')
                
        rmd_bad_evt.close()
        bad_evt.close()


    ############################################## 
    #*-*-*-*-*- Plot path
    #
    # Plot path between event-station pairs
    # 
    
    def plot_path(self):

        nordic_inp = os.path.join('tools','select.out')

        evt_dic = Read_Nordic(nordic_inp)

        vel_dic,sta_dic = self.read_VelSta(inp_vel=os.path.join('velinp',self.velmod_name),inp_sta=self.input_dic['Stationfile'])

        self.sta_list = {}

        stnm  = sta_dic.keys()
        stlat = [sta_dic[i]['lat'] for i in sta_dic]                # station latitude
        stlon = [sta_dic[i]['lon'] for i in sta_dic]                # station longitude


        fig = plt.figure()
        ax  = plt.subplot(1,1,1)

        evrms = []
        evlon = []
        evlat = []
            
        for evt in evt_dic:

            elat  = evt_dic[evt]['HEADER']['L1']['Lat']
            elon  = evt_dic[evt]['HEADER']['L1']['Lon']
            evrms.append(evt_dic[evt]['HEADER']['L1']['RMS'])
            evlon.append(elon)
            evlat.append(elat)

            for sta in evt_dic[evt]['PHASE']:
                
                _ = sta.ljust(4)

                try:

                    slat = sta_dic[_]['lat']
                    slon = sta_dic[_]['lon']
                    ax.plot([elon,slon],[elat,slat],color='grey',linewidth=.7,zorder=1)

                except KeyError:

                    pass
                
        # plot station-event distribution, ray-path

        # plot some features

        self.m        = Basemap(llcrnrlon=self.lon_min, llcrnrlat=self.lat_min,urcrnrlon=self.lon_max, urcrnrlat=self.lat_max,
                                fix_aspect=False)
        self.tehran   = loadtxt(path.join('gdb','db','IRN_PROV','tehran.poly'))
        self.x,self.y = self.tehran[:,0], self.tehran[:,1]
        self.m.plot(self.x,self.y, marker=None,color='r',linewidth=2.0,zorder=1)
        ax.text(51.3890, 35.6892, '$Tehran$')
        ax.plot(51.3700, 35.6892,marker='s', markersize=5, color='k')

        #____ Faults, Borders

        self.m.readshapefile(path.join('gdb','db','IRN_PROV','IRN_adm0'),'iran')
        self.m.readshapefile(path.join('gdb','db','IRN_PROV','IRN_adm1'),'prov', color='black')

        #____ Ticks

        parallels = linspace(self.lat_min, self.lat_max, 5)
        meridians = linspace(self.lon_min, self.lon_max, 5)

        #____ Labels = [left,right,top,bottom]

        parallels = linspace(self.lat_min,self.lat_max,5)
        self.m.drawparallels(parallels,labels=[True,False,False,False],fmt='%4.1f',linewidth=0)
        meridians = linspace(self.lon_min,self.lon_max,5)
        self.m.drawmeridians(meridians,labels=[False,False,False,True],fmt='%4.1f',linewidth=0)

        #____ Scale & North sing

        self.scale_len = self.sc_len # km
        self.sc_lon_st = min(ax.get_xlim())+0.03*diff(ax.get_xlim())[0]
        self.sc_lon_ed = self.sc_lon_st + k2d(self.scale_len)
        self.sc_lat_st = min(ax.get_ylim())+0.03*diff(ax.get_ylim())[0]
        self.sc_lat_ed = self.sc_lat_st
        self.txt_x     = mean([self.sc_lon_st, self.sc_lon_ed])
        self.txt_y     = self.sc_lat_ed + abs(diff([self.sc_lat_ed,min(ax.get_ylim())]))

        ax.plot([self.sc_lon_st,self.sc_lon_ed],
                [self.sc_lat_st,self.sc_lat_ed],
                color='k',linewidth=3.0)
        ax.plot(self.sc_lon_st,self.sc_lat_ed,
                color='k',linewidth=3.0, marker='>',ms=10)
        ax.plot(self.sc_lon_ed,self.sc_lat_ed,
                color='k',linewidth=3.0, marker='<',ms=10)
        ax.text(self.txt_x,self.txt_y,str(self.scale_len)+' km',
                horizontalalignment='center')

        ax.annotate('N', xy=(.030, .07), xycoords='axes fraction', fontsize=16,  ha='center', va='bottom')
        ax.annotate('^', xy=(.031, .09), xycoords='axes fraction', fontsize=22,  ha='center', va='bottom')

        c   = evrms
        im  = ax.scatter(evlon,evlat,c=c,marker='.',s=100,edgecolor='black',cmap=plt.cm.cool,zorder=4,vmin=0.0, vmax=self.max_rms)
        ax.scatter(stlon,stlat, s=200, marker='^', edgecolor='red', linewidth=1.5, facecolor='y',zorder=5,alpha=0.9)
                
        ax.grid(True)
        ax.locator_params(axis='x',nbins=5)
        ax.locator_params(axis='y',nbins=5)
        divider      = make_axes_locatable(ax)
        cax          = divider.append_axes("right", size="2%", pad=0.05)
        cb           = plt.colorbar(im, ax=ax, cax=cax, label='RMS [second]')
        tick_locator = ticker.MaxNLocator(nbins=6)
        cb.locator   = tick_locator
        cb.update_ticks()
        ax.set_xlim(self.lon_min,self.lon_max)
        ax.set_ylim(self.lat_min,self.lat_max)
        ax.set_title('$Path Coverage$')

        plt.savefig(os.path.join('figs',self.project_nm,'path_'+self.project_nm))
        plt.close()


    ############################################## 
    #*-*-*-*-*- Plot final results
    #
    # Make final figures using final location.
    
    def plot_final(self):

        print "+++ Start to plot final figures, Please wait ..."

        def cpt2seg(file_name, sym=False, discrete=False):

            dic    = {}
            f      = open(file_name, 'r')
            rgb    = loadtxt(f)
            rgb    = rgb/255.
            s      = shape(rgb)
            colors = ['red', 'green', 'blue']
            
            for c in colors:
                
                i = colors.index(c)
                x = rgb[:, i+1]

                if discrete:
                    
                    if sym:
                        
                        dic[c]      = zeros((2*s[0]+1, 3), dtype=float)
                        dic[c][:,0] = linspace(0,1,2*s[0]+1)
                        vec         = concatenate((x ,x[::-1]))
                        
                    else:
                        
                        dic[c]      = zeros((s[0]+1, 3), dtype=float)
                        dic[c][:,0] = linspace(0,1,s[0]+1)
                        vec         = x
                        
                    dic[c][1:, 1] = vec
                    dic[c][:-1,2] = vec
                       
                else:
                    
                    if sym:
                        
                        dic[c]      = zeros((2*s[0], 3), dtype=float)
                        dic[c][:,0] = linspace(0,1,2*s[0])
                        vec         = concatenate((x ,x[::-1]))
                        
                    else:
                        
                        dic[c]      = zeros((s[0], 3), dtype=float)
                        dic[c][:,0] = linspace(0,1,s[0])
                        vec         = x
                        
                    dic[c][:, 1] = vec
                    dic[c][:, 2] = vec
           
            return dic


        # read hyp.out file into a dictionary

        def eigsorted(cov):

            vals, vecs = linalg.eigh(cov)
            order      = vals.argsort()[::-1]
            
            return vals[order], vecs[:,order]

        def extract_hyp(hyp_path, region):

            hyp_dic = {}

            nordic_dic = Read_Nordic(hyp_path)

            for evt in sorted(nordic_dic.keys()):

                if nordic_dic[evt]['HEADER']['L1'] and nordic_dic[evt]['HEADER']['LE'] and nordic_dic[evt]['PHASE']:

                    lat = nordic_dic[evt]['HEADER']['L1']['Lat']
                    lon = nordic_dic[evt]['HEADER']['L1']['Lon']

                    if region[2]<=lat<region[3] and region[0]<=lon<region[1]:

                        y    = nordic_dic[evt]['HEADER']['L1']['Year']
                        m    = nordic_dic[evt]['HEADER']['L1']['Month']
                        d    = float(nordic_dic[evt]['HEADER']['L1']['Day'])
                        H    = float(nordic_dic[evt]['HEADER']['L1']['Hour'])
                        M    = float(nordic_dic[evt]['HEADER']['L1']['Min'])
                        S    = float(nordic_dic[evt]['HEADER']['L1']['Sec'])

                        if S >= 60:

                            S = S - 60
                            M+=1

                        if M >= 60:

                            M = M - 60
                            H+=1

                        if H >= 24:

                            H = H - 24
                            d+=1

                        ot   = '%4d%02d%02d %02d%02d%04.1f'%(y,m,d,H,M,S)
                        ot   = dt.strptime(ot,'%Y%m%d %H%M%S.%f')
                        hyp_dic[ot] = {} 
                        hyp_dic[ot]['evlat'] = nordic_dic[evt]['HEADER']['L1']['Lat']
                        hyp_dic[ot]['evlon'] = nordic_dic[evt]['HEADER']['L1']['Lon']
                        hyp_dic[ot]['evdep'] = nordic_dic[evt]['HEADER']['L1']['Dep']
                        hyp_dic[ot]['evmag'] = nordic_dic[evt]['HEADER']['L1']['Mag1']
                        hyp_dic[ot]['evno']  = nordic_dic[evt]['HEADER']['L1']['NumS']
                        hyp_dic[ot]['evdmin'] = []

                        for k in nordic_dic[evt]['PHASE']:

                            for Pphase in nordic_dic[evt]['PHASE'][k]['P']:

                                hyp_dic[ot]['evdmin'].append(nordic_dic[evt]['PHASE'][k]['P'][Pphase]['DIS'])

                        hyp_dic[ot]['evdmin'] = min(hyp_dic[ot]['evdmin'])        
                        hyp_dic[ot]['evrms']  = nordic_dic[evt]['HEADER']['L1']['RMS']
                        
                        if nordic_dic[evt]['HEADER']['LE']['LAE'] and nordic_dic[evt]['HEADER']['LE']['LOE']:

                            hyp_dic[ot]['evgap'] = nordic_dic[evt]['HEADER']['LE']['GAP']
                            hyp_dic[ot]['everh'] = sqrt(nordic_dic[evt]['HEADER']['LE']['LAE']**2 + nordic_dic[evt]['HEADER']['LE']['LOE']**2)
                            hyp_dic[ot]['everx'] = nordic_dic[evt]['HEADER']['LE']['LOE']
                            hyp_dic[ot]['every'] = nordic_dic[evt]['HEADER']['LE']['LAE']
                            hyp_dic[ot]['everz'] = nordic_dic[evt]['HEADER']['LE']['DEE']
                            hyp_dic[ot]['evcxy'] = nordic_dic[evt]['HEADER']['LE']['CXY']
                            cov = array([[hyp_dic[ot]['everx']**2, hyp_dic[ot]['evcxy']], [hyp_dic[ot]['evcxy'], hyp_dic[ot]['every']**2]])

                            vals, vecs            = eigsorted(cov)
                            theta                 = degrees(arctan2(*vecs[:,0][::-1]))
                            width, height         = 2 * self.nstd * sqrt(vals)
                            ellip                 = [(hyp_dic[ot]['evlon'],hyp_dic[ot]['evlat']),k2d(width), k2d(height), theta]
                            hyp_dic[ot]['elipse'] = ellip

            return hyp_dic


        # extract topography from a big SRTM file

        def tif2npz(fid, region, res_x, res_y):

            system("gdalwarp -tr "+res_x+" "+res_y+
                   " -te "+
                   str(region[0])+" "+
                   str(region[2])+" "+
                   str(region[1])+" "+
                   str(region[3])+" "+
                   fid+" cliped.tif > /dev/null")

            ds  = gdal.Open('cliped.tif')
            z   = ds.ReadAsArray()[::-1,:]
            lnv = float(self.bat_dep) # replace large negative values in z
            
            if z.min() <= -1000.0:

                bad_val = z.min()
                
                for i in z:

                    i[i==bad_val] = lnv
                
            savez('cliped.npz',elevation=array(z))
            remove('cliped.tif')


        #___________________ START
            

        #_________ Extract topography

        if self.topo_flag == 'True':

            tif2npz(fid=self.srtm_tif, region=[self.lon_min, self.lon_max, self.lat_min, self.lat_max],
                    res_x=self.res_x, res_y=self.res_y)

            with load("cliped.npz") as dem:
                
                self.elev = dem['elevation']

            remove("cliped.npz")

        #_________ Extract events

        self.iniloc_file = os.path.join('tools','select.out')
        self.finloc_file = os.path.join('tmp',self.project_nm,'report','final','hyp.out')

        self.ini = extract_hyp(hyp_path=self.iniloc_file, region=[self.lon_min, self.lon_max,
                                                             self.lat_min, self.lat_max])
        self.fin = extract_hyp(hyp_path=self.finloc_file, region=[self.lon_min, self.lon_max,
                                                             self.lat_min, self.lat_max])
        self.dbc = 0

        for db in [self.ini,self.fin,[]]:

            self.events = []

            for _ in db:

                self.events.append([db[_]['evlon'],db[_]['evlat'],db[_]['evdep'],db[_]['evmag']])

            if not len(db):

                self.events = []

                for d in self.ini:

                    for e in self.fin:

                        self.tdiff = d-e
                        self.tdiff = abs(self.tdiff.total_seconds())

                        if self.tdiff <= self.max_tdf:

                            self.events.append([self.ini[d]['evlon'],self.ini[d]['evlat'],
                                                self.fin[e]['evlon'],self.fin[e]['evlat']])               

            self.events = array(self.events)
            self.color  = self.events[:,2]    # color --> dep 

            #_________ Plot

            fig = plt.figure()
            ax  = plt.subplot(1,1,1)
            s   = 10

            #____ Main frame

            self.cnt = [self.lon_min+(self.lon_max-self.lon_min)/2.0,
                        self.lat_min+(self.lat_max-self.lat_min)/2.0]

            self.m = Basemap(llcrnrlon=self.lon_min,
                             llcrnrlat=self.lat_min,
                             urcrnrlon=self.lon_max,
                             urcrnrlat=self.lat_max,
                             fix_aspect=False)

            if self.cptflag == 'G':

                self.cmap = cpt2seg(self.gmtcpt)
                self.cmap = LinearSegmentedColormap('palette', self.cmap, 256)

            elif self.cptflag == 'P':

                self.cmap = self.pltcpt

            if self.topo_flag == 'True':

                self.ve   = 0.05
                self.ls   = LightSource(azdeg=135, altdeg=45)
                self.rgb  = self.ls.shade(self.elev, cmap=self.cmap, vert_exag=self.ve, blend_mode='hsv')
                self.m.imshow(self.rgb)

            else:

                self.m.drawmapboundary(fill_color='aqua')
                self.m.fillcontinents(color='#ddaa66',lake_color='aqua')
                    
            #____ Tehran border, symbol and label

            self.tehran      = loadtxt(path.join('gdb','db','IRN_PROV','tehran.poly'))
            self.x, self.y   = self.tehran[:,0], self.tehran[:,1]
            
            self.m.plot(self.x,self.y, marker=None,color='r',linewidth=2.0,zorder=1)
            
            ax.text(51.3890, 35.6892, '$Tehran$')
            ax.plot(51.3700, 35.6892,marker='s', markersize=5, color='k')

            #____ Events

            if self.dbc < 2:

                if self.mag_flag == 'True':
                    
                    for i in range(6):
                        
                        if len(self.events[(self.events[:,3]>i) & (self.events[:,3]<i+1)]):

                            self.tmp = self.events[(self.events[:,3]>i) & (self.events[:,3]<i+1)]
                            self.clr = self.color[(self.events[:,3]>i) & (self.events[:,3]<i+1)]
                            self.x,self.y = self.tmp[:,0],self.tmp[:,1]

                            sc = ax.scatter(self.x,self.y, s=s, c=self.clr,
                                            vmin=self.dep_min, vmax=self.dep_max,
                                           alpha=1.0, cmap=plt.cm.bwr, edgecolor='w',zorder=20)

                        elif len(self.events[(self.events[:,3]>5)]):

                            self.tmp = self.events[(self.events[:,3]>5)]
                            self.clr = self.color[(self.events[:,3]>5)]
                            self.x,self.y = self.tmp[:,0],self.tmp[:,1]

                            im = ax.scatter(self.x,self.y, s=s, c=self.clr,
                                            vmin=self.dep_min, vmax=self.dep_max,
                                           alpha=1.0, cmap=plt.cm.bwr, edgecolor='w',zorder=20)

                        s+=50
                else:

                    self.x,self.y = self.events[:,0],self.events[:,1]

                    im = ax.scatter(self.x,self.y, s=75, c='red',
                                    vmin=self.dep_min, vmax=self.dep_max,
                                   alpha=1.0, cmap=plt.cm.bwr, edgecolor='w',zorder=20)
                
                for _ in db:

                    pos    = self.m(db[_]['elipse'][0][0],db[_]['elipse'][0][1])
                    width  = db[_]['elipse'][1]
                    height = db[_]['elipse'][2]
                    theta  = db[_]['elipse'][3]

                    ellip  = Ellipse(xy=pos, width=width, height=height,angle=theta,
                                     facecolor='none', edgecolor='b', linewidth=1.5,zorder=21)
                    ax.add_artist(ellip)

            else:

                for x0,y0,x1,y1 in self.events:

                    ax.arrow(x0,y0,diff([x0,x1])[0],diff([y0,y1])[0],
                             head_length=0.01,head_width=0.01,width=0.002,fc='b', ec='b')
                  

            #____ Faults

            if self.plt_flt == 'True':

                self.m.readshapefile(path.join('gdb','db','IRN_FLT','fault_iran'),'fault',linewidth=1.5)
                self.m.readshapefile(path.join('gdb','db','IRN_PROV','IRN_adm0'),'iran')
                self.m.readshapefile(path.join('gdb','db','IRN_PROV','IRN_adm1'),'prov', color='grey')

            #____ Ticks

            parallels = linspace(self.lat_min, self.lat_max, 5)
            meridians = linspace(self.lon_min, self.lon_max, 5)

            #____ Labels = [left,right,top,bottom]

            parallels = linspace(self.lat_min,self.lat_max,5)
            self.m.drawparallels(parallels,labels=[True,False,False,False],fmt='%4.1f',linewidth=0)
            self.meridians = linspace(self.lon_min,self.lon_max,5)
            self.m.drawmeridians(meridians,labels=[False,False,False,True],fmt='%4.1f',linewidth=0)

            #____ Scale & North sing

            self.scale_len = self.sc_len # km
            self.sc_lon_st = min(ax.get_xlim())+0.03*diff(ax.get_xlim())[0]
            self.sc_lon_ed = self.sc_lon_st + k2d(self.scale_len)
            self.sc_lat_st = min(ax.get_ylim())+0.03*diff(ax.get_ylim())[0]
            self.sc_lat_ed = self.sc_lat_st
            self.txt_x     = mean([self.sc_lon_st, self.sc_lon_ed])
            self.txt_y     = self.sc_lat_ed + abs(diff([self.sc_lat_ed,min(ax.get_ylim())]))

            ax.plot([self.sc_lon_st,self.sc_lon_ed],
                    [self.sc_lat_st,self.sc_lat_ed],
                    color='k',linewidth=3.0)
            ax.plot(self.sc_lon_st,self.sc_lat_ed,
                    color='k',linewidth=3.0, marker='>',ms=10)
            ax.plot(self.sc_lon_ed,self.sc_lat_ed,
                    color='k',linewidth=3.0, marker='<',ms=10)
            ax.text(self.txt_x,self.txt_y,str(self.scale_len)+' km',
                    horizontalalignment='center')

            ax.annotate('N', xy=(.030, .07), xycoords='axes fraction', fontsize=16,  ha='center', va='bottom')
            ax.annotate('^', xy=(.031, .09), xycoords='axes fraction', fontsize=22,  ha='center', va='bottom')


            #____ Save Figs

            if self.dbc == 0:

                n = 'ini_'+self.project_nm

            elif self.dbc == 1:

                n = 'fin_'+self.project_nm

            else:

                n = 'diff_'+self.project_nm

            ax.set_aspect('equal')

            plt.savefig(path.join('figs',self.project_nm,'topo_'+n+'.png'),dpi=300)
            plt.close()

            self.dbc+=1

        #____ Final plots
            
        fig = plt.figure()

        # plot results, station correction

        self.sta_list = {}

        with open(os.path.join('tmp',self.project_nm,'report','final','STATION0.HYP')) as f:

            for l in f:

                if (l[13:14] == 'N' or l[13:14] == 'S')  and (l[22:23] == 'E' or l[22:23] == 'W'):

                    n   = l[2:6].strip()
                    lat = ll.Latitude(degree=float(l[6:8]), minute=float(l[8:13])).decimal_degree
                    lon = ll.Longitude(degree=float(l[14:17]), minute=float(l[17:22])).decimal_degree
                    
                    self.sta_list[n] = {'lon':lon,'lat':lat}

        self.st_cod = []
        self.st_lat = []
        self.st_lon = []
        self.st_ptc = []
        self.st_stc = []
        self.ev_lat = []
        self.ev_lon = []
        self.ev_rms = []

        for s in self.sta_cor:

            n = s.strip()
            self.st_cod.append(n)
            self.st_lat.append(float(self.sta_list[n]['lat']))
            self.st_lon.append(float(self.sta_list[n]['lon']))
            
            if self.sta_cor[s]['p'].strip():

                self.st_ptc.append(float(self.sta_cor[s]['p']))

            if self.sta_cor[s]['s'].strip():

                self.st_stc.append(float(self.sta_cor[s]['s']))

        self.st_cod = array(self.st_cod)
        self.st_lat = array(self.st_lat)
        self.st_lon = array(self.st_lon)
        self.st_ptc = array(self.st_ptc)
        self.st_stc = array(self.st_stc)

        for _ in self.fin:

            self.ev_lon.append(self.fin[_]['evlon'])
            self.ev_lat.append(self.fin[_]['evlat'])
            self.ev_rms.append(self.fin[_]['evrms'])

        
        self.ev_lon = array(self.ev_lon)
        self.ev_lat = array(self.ev_lat)
        self.ev_rms = array(self.ev_rms)
        

                                  
        ax1 = plt.axes([0.05, 0.45, .45, .5])

        # plot some features

        self.m = Basemap(llcrnrlon=self.lon_min, llcrnrlat=self.lat_min,
                    urcrnrlon=self.lon_max, urcrnrlat=self.lat_max,fix_aspect=False)
        self.m.plot(self.x,self.y, marker=None,color='r',linewidth=2.0,zorder=1)
        ax1.text(51.3890, 35.6892, '$Tehran$')
        ax1.plot(51.3700, 35.6892,marker='s', markersize=5, color='k')


        #____ Faults, Borders

        self.m.readshapefile(path.join('gdb','db','IRN_PROV','IRN_adm0'),'iran')
        self.m.readshapefile(path.join('gdb','db','IRN_PROV','IRN_adm1'),'prov', color='grey')

        #____ Ticks

        parallels = linspace(self.lat_min, self.lat_max, 5)
        meridians = linspace(self.lon_min, self.lon_max, 5)

        #____ Labels = [left,right,top,bottom]

        parallels = linspace(self.lat_min,self.lat_max,5)
        self.m.drawparallels(parallels,labels=[True,False,False,False],fmt='%4.1f',linewidth=0)
        meridians = linspace(self.lon_min,self.lon_max,5)
        self.m.drawmeridians(meridians,labels=[False,False,False,True],fmt='%4.1f',linewidth=0)

        #____ Scale & North sing

        self.scale_len = self.sc_len # km
        self.sc_lon_st = min(ax1.get_xlim())+0.03*diff(ax1.get_xlim())[0]
        self.sc_lon_ed = self.sc_lon_st + k2d(self.scale_len)
        self.sc_lat_st = min(ax1.get_ylim())+0.03*diff(ax1.get_ylim())[0]
        self.sc_lat_ed = self.sc_lat_st
        self.txt_x     = mean([self.sc_lon_st, self.sc_lon_ed])
        self.txt_y     = self.sc_lat_ed + abs(diff([self.sc_lat_ed,min(ax1.get_ylim())]))

        ax1.plot([self.sc_lon_st,self.sc_lon_ed],
                [self.sc_lat_st,self.sc_lat_ed],
                color='k',linewidth=3.0)
        ax1.plot(self.sc_lon_st,self.sc_lat_ed,
                color='k',linewidth=3.0, marker='>',ms=10)
        ax1.plot(self.sc_lon_ed,self.sc_lat_ed,
                color='k',linewidth=3.0, marker='<',ms=10)
        ax1.text(self.txt_x,self.txt_y,str(self.scale_len)+' km',
                horizontalalignment='center')

        ax1.annotate('N', xy=(.030, .07), xycoords='axes fraction', fontsize=16,  ha='center', va='bottom')
        ax1.annotate('^', xy=(.031, .09), xycoords='axes fraction', fontsize=22,  ha='center', va='bottom')


        #_____ Plot events/sta_cor
        
        scl = 1000
        c   = self.ev_rms
        im  = ax1.scatter(self.ev_lon,self.ev_lat,
                          c=c,marker='.',s=100,edgecolor='black',cmap=plt.cm.cool, vmin=0.0, vmax=self.max_rms)
        ax1.scatter(self.st_lon[0],self.st_lat[0],
                    s=0.3*scl,marker='+',color='none',linewidth=1,label='+%.1fsec'%(self.max_rms))
        ax1.scatter(self.st_lon[self.st_ptc>=0],self.st_lat[self.st_ptc>=0],
                    s=abs(self.st_ptc[self.st_ptc>=0]*scl),marker='+',color='blue',linewidth=1)
        ax1.scatter(self.st_lon[0],self.st_lat[0],
                    s=0.3*scl,marker='o',color='none',linewidth=1,label='-%.1fsec'%(self.max_rms))
        ax1.scatter(self.st_lon[self.st_ptc<0] ,self.st_lat[self.st_ptc<0],
                    s=abs(self.st_ptc[self.st_ptc<0]*scl),marker='o',facecolors='none' ,linewidth=1,edgecolor='red')

        for x,y,t in zip (self.st_lon,self.st_lat,self.st_cod):

            if self.lon_min<=x<=self.lon_max and self.lat_min<=y<=self.lat_max:

                ax1.text(x,y,t,fontsize=10)

        plt.grid()
        ax1.set_xlim(self.lon_min,self.lon_max)
        ax1.set_ylim(self.lat_min,self.lat_max)
        ax1.locator_params(axis='x',nbins=5)
        ax1.locator_params(axis='y',nbins=5)
        leg = ax1.legend(loc=4)
        legend = ax1.get_legend()
        legend.legendHandles[0].set_color('red')
        legend.legendHandles[1].set_edgecolor('blue')
        divider      = make_axes_locatable(ax1)
        cax          = divider.append_axes("right", size="2%", pad=0.05)
        cb           = plt.colorbar(im, ax=ax1, cax=cax, label='RMS [second]')
        tick_locator = ticker.MaxNLocator(nbins=6)
        cb.locator   = tick_locator
        cb.update_ticks()


        self.ini_lon = []
        self.ini_lat = []
        self.ini_dep = []
        self.ini_rms = []
        self.ini_erh = []
        self.ini_erz = []
        self.fin_lon = []
        self.fin_lat = []
        self.fin_dep = []
        self.fin_rms = []
        self.fin_erh = []
        self.fin_erz = []

        for i in sorted(self.ini):

            i_lon = self.ini[i]['evlon']
            i_lat = self.ini[i]['evlat'] 

            for f in sorted(self.fin):

                f_lon  = self.fin[f]['evlon']
                f_lat  = self.fin[f]['evlat']
                dx     = f_lon - i_lon
                dy     = f_lat - i_lat
                diff_h = sqrt(dx**2 + dy**2)
                
                diff_ot = i-f
                diff_ot = abs(diff_ot.total_seconds())
                

                if diff_ot <= self.max_tdf and diff_h <= self.max_ldf:

                    self.ini_lon.append(self.ini[i]['evlon'])
                    self.ini_lat.append(self.ini[i]['evlat'])
                    self.ini_dep.append(self.ini[i]['evdep'])
                    self.ini_rms.append(self.ini[i]['evrms'])
                    self.ini_erh.append(self.ini[i]['everh'])
                    self.ini_erz.append(self.ini[i]['everz'])

                    self.fin_lon.append(self.fin[f]['evlon'])
                    self.fin_lat.append(self.fin[f]['evlat'])
                    self.fin_dep.append(self.fin[f]['evdep'])
                    self.fin_rms.append(self.fin[f]['evrms'])
                    self.fin_erh.append(self.fin[f]['everh'])
                    self.fin_erz.append(self.fin[f]['everz'])

                    break


        self.ini_lon = array(self.ini_lon)
        self.ini_lat = array(self.ini_lat)
        self.ini_dep = array(self.ini_dep)
        self.ini_rms = array(self.ini_rms)
        self.ini_erh = array(self.ini_erh)
        self.ini_erz = array(self.ini_erz)
        self.fin_lon = array(self.fin_lon)
        self.fin_lat = array(self.fin_lat)
        self.fin_dep = array(self.fin_dep)
        self.fin_rms = array(self.fin_rms)
        self.fin_erh = array(self.fin_erh)
        self.fin_erz = array(self.fin_erz)

        
        # Plot Latitude Dislocation
        
        ax2 = plt.axes([0.60, 0.70, .37, .25])

        for (y,yy) in zip(self.ini_lat,self.fin_lat):
            if yy-y>0:
                ax2.plot([y,y],[d2k(yy-y),0],'r-o',markevery=self.ini_lat.size)
            else:
                ax2.plot([y,y],[d2k(yy-y),0],'b-o',markevery=self.ini_lat.size)
        ax2.set_xlabel('Latitude [deg]')
        ax2.set_ylabel('Dislocation [km]')
        ax2.set_ylim(-max([abs(_) for _ in ax2.get_ylim()]),max([abs(_) for _ in ax2.get_ylim()]))
        ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        xlim=ax2.get_xlim()
        ax2.plot(xlim,[0,0], linewidth=2, color = 'k')
        ax2.set_xlim(xlim)
        ax2.locator_params(axis='x',nbins=5)
        ax2.locator_params(axis='y',nbins=5)
        plt.grid()

        # Plot Longitude Dislocation
        
        ax3 = plt.axes([0.60, 0.38, .37, .25])
        
        for (x,xx) in zip(self.ini_lon,self.fin_lon):
            
            if xx-x>0:
                ax3.plot([x,x],[d2k(xx-x),0],'r-o',markevery=self.ini_lon.size)
            else:
                ax3.plot([x,x],[d2k(xx-x),0],'b-o',markevery=self.ini_lon.size)
                
        ax3.set_xlabel('Longitude [deg]')
        ax3.set_ylabel('Dislocation [km]')
        ax3.set_ylim(-max([abs(_) for _ in ax3.get_ylim()]),max([abs(_) for _ in ax3.get_ylim()]))
        ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        xlim=ax3.get_xlim()
        ax3.plot(xlim,[0,0], linewidth=2, color = 'k')
        ax3.set_xlim(xlim)
        ax3.locator_params(axis='x',nbins=5)
        ax3.locator_params(axis='y',nbins=5)
        plt.grid()

        # Plot Depth Dislocation
        
        ax4 = plt.axes([0.60, 0.06, .37, .25])
        
        for (x,xx) in zip(self.ini_dep,self.fin_dep):
            
            if xx-x>0:
                ax4.plot([x,x],[xx-x,0],'r-o',markevery=self.ini_dep.size)
            else:
                ax4.plot([x,x],[xx-x,0],'b-o',markevery=self.ini_dep.size)
                
        ax4.set_xlabel('Depth [km]')
        ax4.set_ylabel('Dislocation [km]')
        ax4.set_ylim(-max([abs(_) for _ in ax4.get_ylim()]),max([abs(_) for _ in ax4.get_ylim()]))
        ax4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax4.plot([0,self.dep_max],[0,0], linewidth=2, color = 'k')
        ax4.set_xlim(0,self.dep_max)
        ax4.locator_params(axis='x',nbins=5)
        ax4.locator_params(axis='y',nbins=5)
        plt.grid()
        
        # Plot Horizontal Error
        
        ax5 = plt.axes([0.05, 0.06, .14, .3])
        bns = arange(0, self.max_her + 1, self.bin_her)
        ax5.hist(self.ini_erh,color='r',label='Initial',alpha=.65,bins=bns)
        ax5.hist(self.fin_erh,color='b',label='Final',alpha=.60,bins=bns)
        ax5.set_xlabel('H_Error [km]')
        ax5.set_ylabel('Event [#]')
        ax5.set_xlim(0,self.max_her)
        ax5.locator_params(axis='x',nbins=5)
        ax5.locator_params(axis='y',nbins=5)
        plt.grid()
        ax5.legend(loc=1,fontsize=plt.rcParams['font.size'])

        # Plot Depth Error
        
        ax6 = plt.axes([0.23, 0.06, .14, .3])
        bns = arange(0, self.max_der + 1, self.bin_der)
        ax6.hist(self.ini_erz,color='r',label='Initial',alpha=.65,bins=bns)
        ax6.hist(self.fin_erz,color='b',label='Final',alpha=.60,bins=bns)
        ax6.set_xlabel('Z_Error [km]')
        ax6.set_xlim(0,self.max_der)
        ax6.locator_params(axis='x',nbins=5)
        ax6.locator_params(axis='y',nbins=5)
        plt.grid()
        ax6.legend(loc=1,fontsize=plt.rcParams['font.size'])

        # Plot RMS
        
        ax7 = plt.axes([0.41, 0.06, .14, .3])
        bns = arange(0, self.max_rms + self.max_rms/20.0, self.max_rms/20.0)
        ax7.hist(self.ini_rms,color='r',label='Initial',alpha=.65,bins=bns)
        ax7.hist(self.fin_rms,color='b',label='Final',alpha=.60,bins=bns)
        ax7.set_xlabel('RMS [sec]')
        ax7.set_xlim(0,self.max_rms)
        ax7.locator_params(axis='x',nbins=5)
        ax7.locator_params(axis='y',nbins=5)
        plt.grid()
        ax7.legend(loc=1,fontsize=plt.rcParams['font.size'])

        plt.savefig(os.path.join('figs',self.project_nm,'final_'+self.project_nm))
        plt.close()

        remove('velest.cmn')
        
        print "+++ Finito!"

#___________________ RUN

import warnings
warnings.filterwarnings("ignore")

start = main()
start.extract_topo()
start.clean()
start.write_par()
start.plot_data()
start.plot_path()
start.mk_syn_model()
start.run_velest()
start.analyse_velest_res()
start.mk_report()
start.plot_final()
