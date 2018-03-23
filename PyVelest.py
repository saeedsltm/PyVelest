#___________________IMPORT FROM NUMPY, SCIPY
from numpy import loadtxt, diff, ones, random, array, hstack, mean, pi, where, median, degrees, ones_like
from numpy import arange, sqrt, linalg, arccos, arctan2, rad2deg, zeros_like, abs, linspace, savez, load
from scipy import zeros, shape, concatenate
#___________________IMPORT EXTRA LIBRARIES
##from osgeo import gdal
import tarfile
import operator
import sys
import glob
from os import system, remove, path, getcwd, chdir, mkdir, makedirs, sep
from datetime import datetime as dt
sys.path.append(path.join('tools','utility'))
from PyNordicRW import Read_Nordic
from PyNordic2CNV import CNV2Nordic, ReadCNV
from collections import OrderedDict
from shutil import copy, rmtree, move
from LatLon import lat_lon as ll
#___________________IMPORT FROM MATPLOTLIB
import mpl_toolkits
mpl_toolkits.__path__.append('/usr/lib/python2.7/dist-packages/mpl_toolkits/')
from initial_mpl import init_plotting_isi
from matplotlib.colors import LightSource
from matplotlib.cbook import get_sample_data
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import ticker, gridspec
import matplotlib.pyplot as plt

"""
ChangeLog:

2018-03-23: Fixed bug for velest backup result (added line 746)

"""

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
        self.input_dic  = self.read_par()

    #*-*-*-*-*- Cheeking if topography data exists
    #
    # Extract tar.gz file may take times

    def extract_topo(self):

        if self.topo_flag == 'True':

            here = getcwd()
            chdir(path.join('gdb','db','srtm'))

            if not path.exists('IRN.tif'):

                print "+++ It seems, You're running 'PyVelest' for the first time on your system!"
                print "+++ You'll never see this message next time."
                print "\n--> Extracting topography data , Please wait ...\n"
                
                tar = tarfile.open('IRN.tif.tar.gz', "r:gz")
                tar.extractall()
                tar.close()

            chdir(here)
    
    #*-*-*-*-*- Remove old files and do empy all folders
    #
    # Remove files from velout.
    # Remove tmp folder if exist.

    def clean(self):

        if path.exists(path.join('tmp',self.project_nm)):

            rmtree(path.join('tmp',self.project_nm))

        if path.exists(path.join('velout',self.project_nm)):

            rmtree(path.join('velout',self.project_nm))

    def write_cmn(self, input_dic, output='velest.cmn', project_nm='test'):

        output     = open(output,'w')
        now        = dt.now()

        header1 = "CONTROL-FILE FOR PROGRAM  V E L E S T  (%s)"%(now.strftime('%Y-%m-%d'))
        header1 = header1.center(70,'*')
        header2 = project_nm

        # write Headers
        
        output.write(header1+'\n')
        output.write(header2+'\n')

        # write the Main lines

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
        output.write(path.join('velinp',self.velmod_name)+'\n')
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

    def read_par(self, par_file=path.join('par','par.dat')):   

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
        self.plt_stn               = par_dic['PLT_STN']
        self.ot_adj                = float(par_dic['OT_ADJ'])
        self.hr_adj                = float(par_dic['HR_ADJ'])
        self.dp_adj                = float(par_dic['DP_ADJ'])

        # Read velocity models in velinp directory and let user to choose wich one must be used

        inp_syntvel      = path.join('par','syntvel.dat')
        inp              = loadtxt(inp_syntvel,dtype=str,comments='#',delimiter=':')
        self.velmod_name = inp[16][1].strip()

        return self.input_dic
    
    #*-*-*-*-*- write new VELEST parameter to velest.cmn file
    #
    # Write VELEST parameters into VELEST's cnv
    # file [./velinp/velest.cmn]
    
    def write_par(self):

        self.project_nm = self.write_cmn(input_dic=self.input_dic, output='velest.cmn',project_nm=self.project_nm)

        print "+++ Write VELEST control file ..."
                            
    #*-*-*-*-*- read velocity model and station files
    #
    # Read VELEST velocity file [./velinp/model.mod] and
    # station file [./velinp/station.sta] and parse it into
    # a dictionary.

    def read_VelSta(self,inp_vel=path.join('velinp','model.mod'), inp_sta=path.join('velinp','station.sta')):

        inp_vel=path.join('velinp',self.velmod_name)

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

    def mk_syn_model(self,inp_syntvel=path.join('par','syntvel.dat')):

        print "+++ Make synthetic velocity models ..."

        inp       = loadtxt(inp_syntvel,dtype=str,comments='#',delimiter=':')
        
        vel_dic,_ = self.read_VelSta(inp_vel=path.join('velinp',self.velmod_name),inp_sta=self.input_dic['Stationfile'])
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
       
        
        if not path.exists(path.join('tmp',self.project_nm,'synthvel')):

            makedirs(path.join('tmp',self.project_nm,'synthvel'))

        else:

            rmtree(path.join('tmp',self.project_nm,'synthvel'))
            makedirs(path.join('tmp',self.project_nm,'synthvel'))
            
        output      = path.join('tmp',self.project_nm,'synthvel')
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
            
            fid = file(path.join(output,"model_"+str(i+1)+".mod"),'w')
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

        init_plotting_isi(9,9)
        plt.rcParams['axes.labelsize'] = 8
        
        ax = plt.subplot(111)
        [i.set_linewidth(0.6) for i in ax.spines.itervalues()]

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
                
                ax.plot(x, y, 'r-', linewidth=1.0, label='Random')
                
            else:
                
                ax.plot(x, y, 'r-', linewidth=1.0)
                ax.set_ylim(ymin=-self.dep_max,ymax=-self.dep_min)
                ax.set_title("Synthetic Models")
                ax.set_xlabel("Velocity (km/s)")
                ax.set_ylabel("Depth (km)")

        velmod_org = [xx,yy]
        
        ax.plot(xx, yy, 'k-', linewidth=1.0,label='Reference')
        ax.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
        ax.set_xlim(self.vel_min, self.vel_max)
        ax.locator_params(axis='x',nbins=7)
        ax.locator_params(axis='y',nbins=7)
        plt.legend(loc=1, fontsize=6)
        
        ax_ = inset_axes(ax, width="30%",  height="40%", loc=3)
        [i.set_linewidth(0.6) for i in ax_.spines.itervalues()]
        ax_.tick_params(labelsize=6)
        binwidth = 2.0
        data     = self.orig_dep
        bins     = arange(min(data), max(data) + binwidth, binwidth)
        ax_.hist(self.orig_dep, bins=bins, orientation="horizontal", rwidth=0.8, color='grey', edgecolor='r')
        ax_.locator_params(axis='x',nbins=5)
        ax_.locator_params(axis='y',nbins=5)
        ax_.set_xlabel('Number of EQs (#)', fontsize=6)
        ax_.xaxis.set_label_position("top")
        ax_.xaxis.tick_top()
        ax_.set_ylabel('Depth (km)', fontsize=6)
        ax_.yaxis.set_label_position("right")
        ax_.yaxis.tick_right()
        ax_.set_ylim(0,self.dep_max)
        ax_.set_ylim(ax_.get_ylim()[::-1])
        ax_.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
        
        plt.tight_layout()
        plt.savefig(path.join('figs',self.project_nm,"synt_models_"+self.project_nm+'.tiff'))
        plt.close()

    #*-*-*-*-*- Run VELEST program
    #
    # Run VELEST for each model
    # 

    def run_velest(self):

        print "+++ Run VELEST ..."

        if not path.exists(path.join('velout',self.project_nm)):

            makedirs(path.join('velout',self.project_nm))

        if path.exists(path.join('tmp',self.project_nm,'synthvel')):

            self.synthvel_files = glob.glob(path.join('tmp',self.project_nm,'synthvel','*.mod'))
            self.synthvel_files = sorted(self.synthvel_files, key=lambda x: int(x.split(sep)[-1].split('_')[1].split('.')[0]))

            print '   --- %d synthetic model(s) found ...'%(len(self.synthvel_files))

        else:

            self.synthvel_files = []
            self.nmod           = 0

        if self.synthvel_files:

            copy('velest.cmn',path.join('tmp',self.project_nm))
    
            for m,c in zip(self.synthvel_files,range(len(self.synthvel_files))):

                tmp_cmn = open(path.join('tmp',self.project_nm,'tmp.cmn'),'w')
                
                with open(path.join('tmp',self.project_nm,'velest.cmn')) as f:

                    for line in f:

                        if path.join('velinp',self.velmod_name) in line:

                            tmp_cmn.write(line.replace(path.join('velinp',self.velmod_name), m))

                        elif path.join('velout','velest.out') in line:

                            tmp_cmn.write(line.replace(path.join('velout','velest.out'),
                                                       path.join('velout',self.project_nm,'velest'+str(c+1)+'.out')))

                        elif path.join('velout','final_loc.cnv') in line:

                            tmp_cmn.write(line.replace(path.join('velout','final_loc.cnv'),
                                                       path.join('velout',self.project_nm,'final_loc'+str(c+1)+'.cnv')))

                        elif path.join('velout','stations_corr.sta') in line:
                            tmp_cmn.write(line.replace(path.join('velout','stations_corr.sta'),
                                                       path.join('velout',self.project_nm,'stations_corr'+str(c+1)+'.sta')))

                        else:

                            tmp_cmn.write(line)

                tmp_cmn.close()
                
                move(path.join('tmp',self.project_nm,'tmp.cmn'),'velest.cmn')

                print '   --- Run for synthetic model:',c+1

                system('velest > /dev/null')

            remove(path.join('tmp',self.project_nm,'velest.cmn'))

            for _ in [path.join('velout','velest.out'),
                      path.join('velout','final_loc.cnv'),
                      path.join('velout','stations_corr.sta')]:

                if path.exists(_):

                    remove(_)

        else:

            print '   --- No synthetic model was found ...'
            print '   --- Run for initial velocity model ...'
            system('velest > /dev/null')

    #*-*-*-*-*- Analyse VELEST results 
    #
    # Analyse VELEST results and plot 
    # them.
    
    def analyse_velest_res(self):

        print "+++ Analyse VELEST outputs ..."

        inp_syntvel    = path.join('par','syntvel.dat')
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
                
                with open(path.join('velout',self.project_nm,'velest%d.out'%(sm+1))) as f:

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

                        if '~~~ solution is worse --> backup ...' in l and iter_counter <= self.nit:

                            self.rms_res[sm].pop(-1)
                            self.avg_abs_adj[sm].pop(-1)
                            iter_counter-=1
                            self.p_vel_adj_flag = False
                            self.s_vel_adj_flag = False
                            
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

            init_plotting_isi(18,10)
            plt.rcParams['axes.labelsize'] = 7
        
            output     = path.join('figs','invres_'+self.project_nm)
            iter_color = ['k','r','g','b','c','y','m','purple','pink','brown'] # iteration color, maximum 10
                               
            ax1         = plt.subplot(2,3,1)
            [i.set_linewidth(0.6) for i in ax1.spines.itervalues()]
            self.minmod = {} 
            for nm in self.p_good_models:
                x      = range(1,1+self.nit)
                ax1.plot(x, [self.avg_abs_adj[nm][_][0] for _ in range(len(self.avg_abs_adj[nm]))], marker='o', ms=3, lw=.75,color='grey')
                self.minmod[nm] = mean([self.avg_abs_adj[nm][_][0] for _ in range(len(self.avg_abs_adj[nm]))])
            self.ind = sorted(self.minmod.items(), key=operator.itemgetter(1))[0][0]
            self.minmod = [self.avg_abs_adj[self.ind][i][0] for i in range(len(self.avg_abs_adj[self.ind]))]
            ax1.plot(x, self.minmod,marker='o',color='red', ms=3, lw=.75,label='Minimum model,%d'%(self.ind))
            ax1.set_xlabel('Iteration (#)')
            ax1.set_ylabel('Origin Time Adjustment (s)')
            ax1.set_ylim(0.0,self.ot_adj)
            ax1.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
            ax1.locator_params(axis='x',nbins=5)
            ax1.locator_params(axis='y',nbins=5)
            ax1.legend(loc=1, fontsize=5)
            ax1.set_xticks(x)
            ax1.set_xticklabels(x)

            ax2 = plt.subplot(2,3,2)
            [i.set_linewidth(0.6) for i in ax2.spines.itervalues()]
            self.minmod = {}
            for nm in self.p_good_models:
                x      = range(1,1+self.nit)
                ax2.plot(x, [self.avg_abs_adj[nm][_][1] for _ in range(len(self.avg_abs_adj[nm]))], marker='o', ms=3, lw=.75, color='grey')
                self.minmod[nm]=mean([self.avg_abs_adj[nm][_][1] for _ in range(len(self.avg_abs_adj[nm]))])
            self.ind = sorted(self.minmod.items(), key=operator.itemgetter(1))[0][0]
            self.minmod = [self.avg_abs_adj[self.ind][i][1] for i in range(len(self.avg_abs_adj[self.ind]))]
            ax2.plot(x, self.minmod,marker='o', ms=3, lw=.75, color='red',label='Minimum model,%d'%(self.ind))
            ax2.set_xlabel('Iteration (#)')
            ax2.set_ylabel('Longitude Adjustment (km)')
            ax2.set_ylim(0.0,self.hr_adj)
            ax2.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
            ax2.locator_params(axis='x',nbins=5)
            ax2.locator_params(axis='y',nbins=5)
            ax2.legend(loc=1, fontsize=5)
            ax2.set_xticks(x)
            ax2.set_xticklabels(x)
            
            ax3 = plt.subplot(2,3,3)
            [i.set_linewidth(0.6) for i in ax3.spines.itervalues()]
            self.minmod = {}
            for nm in self.p_good_models:
                x      = range(1,1+self.nit)
                ax3.plot(x, [self.avg_abs_adj[nm][_][2] for _ in range(len(self.avg_abs_adj[nm]))], marker='o', ms=3, lw=.75, color='grey')
                self.minmod[nm]=mean([self.avg_abs_adj[nm][_][2] for _ in range(len(self.avg_abs_adj[nm]))])
            self.ind = sorted(self.minmod.items(), key=operator.itemgetter(1))[0][0]
            self.minmod = [self.avg_abs_adj[self.ind][i][2] for i in range(len(self.avg_abs_adj[self.ind]))]
            ax3.plot(x, self.minmod,marker='o', ms=3, lw=.75, color='red',label='Minimum model,%d'%(self.ind))
            ax3.set_xlabel('Iteration (#)')
            ax3.set_ylabel('Latitude Adjustment (km)')
            ax3.set_ylim(0.0,self.hr_adj)
            ax3.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
            ax3.locator_params(axis='x',nbins=5)
            ax3.locator_params(axis='y',nbins=5)
            ax3.legend(loc=1, fontsize=5)
            ax3.set_xticks(x)
            ax3.set_xticklabels(x)
            
            ax4 = plt.subplot(2,3,4)
            [i.set_linewidth(0.6) for i in ax4.spines.itervalues()]
            self.minmod = {}
            for nm in self.p_good_models:
                x      = range(1,1+self.nit)
                ax4.plot(x, [self.avg_abs_adj[nm][_][3] for _ in range(len(self.avg_abs_adj[nm]))], marker='o', ms=3, lw=.75, color='grey')
                self.minmod[nm]=mean([self.avg_abs_adj[nm][_][3] for _ in range(len(self.avg_abs_adj[nm]))])
            self.ind = sorted(self.minmod.items(), key=operator.itemgetter(1))[0][0]
            self.minmod = [self.avg_abs_adj[self.ind][i][3] for i in range(len(self.avg_abs_adj[self.ind]))]
            ax4.plot(x, self.minmod,marker='o', ms=3, lw=.75, color='red',label='Minimum model,%d'%(self.ind))
            ax4.set_xlabel('Iteration (#)')
            ax4.set_ylabel('Depth Adjustment (km)')
            ax4.set_ylim(0.0,self.dp_adj)
            ax4.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
            ax4.locator_params(axis='x',nbins=5)
            ax4.locator_params(axis='y',nbins=5)
            ax4.legend(loc=1, fontsize=5)
            ax4.set_xticks(x)
            ax4.set_xticklabels(x)
            
            ax5 = plt.subplot(2,3,5)
            [i.set_linewidth(0.6) for i in ax5.spines.itervalues()]
            self.minmod = {}
            self.test = {}
            for nm in self.p_good_models:
                x      = range(0,1+self.nit)
                ax5.plot(x, self.rms_res[nm],marker='o', ms=3, lw=.75, color='grey')
                self.minmod[nm]=mean(self.rms_res[nm])
            self.ind = sorted(self.minmod.items(), key=operator.itemgetter(1))[0][0]
            self.minmod = self.rms_res[self.ind]
            ax5.plot(x, self.minmod,marker='o', ms=3, lw=.75, color='red',label='Minimum model,%d'%(self.ind))
            ax5.set_xlabel('Iteration (#)')
            ax5.set_ylabel('RMS-Residual (s)')
            ax5.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
            ax5.locator_params(axis='x',nbins=5)
            ax5.locator_params(axis='y',nbins=5)
            ax5.legend(loc=1, fontsize=5)
            ax5.set_xticks(x)
            ax5.set_xticklabels(x)
            
            ax6 = plt.subplot(2,3,6)
            [i.set_linewidth(0.6) for i in ax6.spines.itervalues()]

            x,y, ymax = [],[],[]

            for it,c,s in zip(self.vel_adj[self.ind],range(len(iter_color)),linspace(3,12,self.nit)[::-1]):
                
                color = iter_color[c]
                                
                for p in it:
                    
                    x.append(p[0])
                    y.append(p[1])
                    
                ymax.append(max([abs(_) for _ in y]))
                it = array(it)
                ax6.plot([it[:,1],zeros_like(it[:,1])],[it[:,0],it[:,0]],
                         markevery=it[:,0].size,lw=s,c=color,label='Iter: '+str(c+1))

            ax6.set_xlabel('Velocity Adjusment (km/s)')
            ax6.set_ylabel('Depth (km)')
            ax6.set_xlim(-max(ymax),max(ymax))
            ax6.locator_params(axis='x',nbins=5)
            ax6.locator_params(axis='y',nbins=5)
            ax6.set_ylim(ax6.get_ylim()[::-1])
            ax6.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            lgnd=plt.legend(by_label.values(), by_label.keys(),fontsize=5)
            for legobj in lgnd.legendHandles:
                legobj.set_linewidth(3.0)
            ax6.axvline(lw=1,color='k',linestyle='--')

            plt.tight_layout()
            plt.savefig(path.join('figs',self.project_nm,'velest_analytics_'+self.project_nm+'.tiff'))
            plt.close()

            # plot results, final velocity, station correction beside final location using minimum model
            
            cnv_dic         = ReadCNV(inpcnv=path.join('velout',self.project_nm,'final_loc%d.cnv'%(self.ind+1)))
            vel_dic,sta_dic = self.read_VelSta(inp_vel=path.join('velinp',self.velmod_name),inp_sta=self.input_dic['Stationfile'])

            print "+++ Plotting data statistics ..."

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

            gs  = gridspec.GridSpec(4, 4)
            fig = plt.figure()

            init_plotting_isi(17,15)
            plt.rcParams['axes.labelsize'] = 7

            ax1 = fig.add_subplot(gs[0:2,0:2])
            [i.set_linewidth(0.6) for i in ax1.spines.itervalues()]
            
            for nm in self.p_good_models:
                ax1.plot(self.svel_p[nm],self.sdep_p[nm],color='red',linewidth=1.0)
            ax1.plot(self.svel_p[0],self.sdep_p[0],color='red',linewidth=1.0,label='Inverted')
            ax1.plot(ovel,odep,color='black',linewidth=1.0,label='Refrence')
            ax1.plot(self.svel_p[self.ind],self.sdep_p[self.ind],color='green',linewidth=1.0,label='Minimum')
            ax1.plot(mean(self.svel_p,axis=0),mean(self.sdep_p,axis=0),color='blue',linewidth=1.0,label='Mean')
            xlim = ax1.get_xlim()
            ylim = ax1.get_ylim()
            xlim = (self.vel_min,self.vel_max)
            ylim = [min(odep),self.dep_max]
            ax1.set_xlim(xlim)
            ax1.set_ylim(ylim)
            ax1.set_ylim(ax1.get_ylim()[::-1])
            ax1.set_xlabel('Velocity (km/s)')
            ax1.set_ylabel('Depth (km)')
            ax1.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
            ax1.legend(numpoints=1, loc=1, fontsize=5)
            ax1.locator_params(axis='x',nbins=8)
            ax1.locator_params(axis='y',nbins=5)

            ax1_ = inset_axes(ax1, width="30%",  height="40%", loc=3)
            [i.set_linewidth(0.6) for i in ax1_.spines.itervalues()]
            ax1_.tick_params(labelsize=6)
        
            binwidth = 2.0
            data     = self.orig_dep
            bins     = arange(min(data), max(data) + binwidth, binwidth)
            ax1_.hist(self.orig_dep, bins=bins, orientation="horizontal", rwidth=0.8, color='grey', edgecolor='r')
            ax1_.locator_params(axis='x',nbins=5)
            ax1_.locator_params(axis='y',nbins=5)
            ax1_.set_xlabel('Number of EQs (#)', fontsize=5)
            ax1_.xaxis.set_label_position("top")
            ax1_.xaxis.tick_top()
            ax1_.set_ylabel('Depth (km)', fontsize=5)
            ax1_.yaxis.set_label_position("right")
            ax1_.yaxis.tick_right()
            ax1_.set_ylim(0,self.dep_max)
            ax1_.set_ylim(ax1_.get_ylim()[::-1])
            ax1_.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
            
            # plot results, station correction

            stnm  = []
            stlat = []
            stlon = []
            stcor = []

            with open(path.join('velout',self.project_nm,'stations_corr%d.sta'%(self.ind+1))) as f:

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
            [i.set_linewidth(0.6) for i in ax2.spines.itervalues()]

            # plot some features

            self.m = Basemap(llcrnrlon=self.lon_min, llcrnrlat=self.lat_min,
                        urcrnrlon=self.lon_max, urcrnrlat=self.lat_max,fix_aspect=False)
            self.tehran = loadtxt(path.join('gdb','db','IRN_PROV','tehran.poly'))
            self.x,self.y = self.tehran[:,0], self.tehran[:,1]
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
            self.m.drawparallels(parallels,labels=[True,False,False,False],fmt='%4.1f',linewidth=0, fontsize=6)
            meridians = linspace(self.lon_min,self.lon_max,5)
            self.m.drawmeridians(meridians,labels=[False,False,False,True],fmt='%4.1f',linewidth=0, fontsize=6)

            #____ Scale & North sing

            self.scale_len = self.sc_len # km
            self.sc_lon_st = min(ax2.get_xlim())+0.04*diff(ax2.get_xlim())[0]
            self.sc_lon_ed = self.sc_lon_st + k2d(self.scale_len)
            self.sc_lat_st = min(ax2.get_ylim())+0.04*diff(ax2.get_ylim())[0]
            self.sc_lat_ed = self.sc_lat_st
            self.txt_x     = mean([self.sc_lon_st, self.sc_lon_ed])
            self.txt_y     = self.sc_lat_ed + abs(diff([self.sc_lat_ed,min(ax2.get_ylim())]))

            ax2.plot([self.sc_lon_st,self.sc_lon_ed],
                    [self.sc_lat_st,self.sc_lat_ed],
                    color='k',linewidth=1.0)
            ax2.plot(self.sc_lon_st,self.sc_lat_ed,
                    color='k',linewidth=1.0, marker='>',ms=4)
            ax2.plot(self.sc_lon_ed,self.sc_lat_ed,
                    color='k',linewidth=1.0, marker='<',ms=4)
            ax2.text(self.txt_x,self.txt_y,str(self.scale_len)+' km',
                    horizontalalignment='center')

            ax2.annotate('N', xy=(.030, .07), xycoords='axes fraction', fontsize=8,  ha='center', va='bottom')
            ax2.annotate('^', xy=(.031, .09), xycoords='axes fraction', fontsize=12,  ha='center', va='bottom')

            dy = ax2.get_ylim()
            dy = dy[1] - dy[0]
            dy = dy
            fs = 0.8*plt.rcParams['font.size']
            scl = 250
            c   = evrms
            im  = ax2.scatter(evlon,evlat,c=c,marker='.',s=50,edgecolor='black',cmap=plt.cm.jet,vmin=0.0, vmax=self.max_rms)
            ax2.scatter(stlon[0],stlat[0],s=0.5*scl,marker='+',color='none',linewidth=1,label='+%.1f'%(self.max_rms))
            ax2.scatter(stlon[stcor>=0],stlat[stcor>=0],s=abs(stcor[stcor>=0]*scl),marker='+',color='blue',linewidth=1)
            ax2.scatter(stlon[0],stlat[0],s=0.5*scl,marker='o',color='none',linewidth=1,label='-%.1f'%(self.max_rms))
            ax2.scatter(stlon[stcor<0] ,stlat[stcor<0], s=abs(stcor[stcor<0]*scl),marker='o',facecolors='none' ,linewidth=1,edgecolor='red')

            for x,y,t in zip (stlon,stlat,stnm):

                if self.lon_min<=x<=self.lon_max and self.lat_min<=y<=self.lat_max:

                    ax2.text(x,y,t,fontsize=fs)

            ax2.set_xlim(self.lon_min,self.lon_max)
            ax2.set_ylim(self.lat_min,self.lat_max)
            ax2.locator_params(axis='x',nbins=5)
            ax2.locator_params(axis='y',nbins=5)
            leg = ax2.legend(loc=4)
            legend = ax2.get_legend()
            legend.legendHandles[0].set_color('blue')
            legend.legendHandles[1].set_edgecolor('red')
            divider      = make_axes_locatable(ax2)
            cax          = divider.append_axes("right", size="4%", pad=0.05)
            cb           = plt.colorbar(im, ax=ax2, cax=cax, label='RMS (s)')
            tick_locator = ticker.MaxNLocator(nbins=6)
            cb.locator   = tick_locator
            cb.update_ticks()
            cb.outline.set_linewidth(.75)


            # plot compare figure between original and minimum model

            self.orig_cnv = ReadCNV(inpcnv=self.input_dic['EQ_datafile'])
            self.mmod_cnv = ReadCNV(inpcnv=path.join('velout',self.project_nm,'final_loc%d.cnv'%(self.ind+1)))

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
            [i.set_linewidth(0.6) for i in ax3.spines.itervalues()]
            
            for (x,xx) in zip(olon,mlon):
                
                if xx-x>0:

                    ax3.plot([x,x],[d2k(xx-x),0],'r-o',markevery=olon.size, ms=3, lw=.5, mew=.1, mec='k')

                else:

                    ax3.plot([x,x],[d2k(xx-x),0],'b-o',markevery=olon.size, ms=3, lw=.5, mew=.1, mec='k')
                    
            ax3.set_xlabel('Longitude')
            ax3.set_ylabel('Dislocation (km)')
            ax3.set_ylim(-max([abs(_) for _ in ax3.get_ylim()]),max([abs(_) for _ in ax3.get_ylim()]))
            ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            xlim=ax3.get_xlim()
            ax3.plot(xlim,[0,0], linewidth=1, color='k')
            ax3.set_xlim(xlim)
            ax3.locator_params(axis='x',nbins=5)
            ax3.locator_params(axis='y',nbins=5)
            ax3.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)

            ax4 = fig.add_subplot(gs[2:,1])
            [i.set_linewidth(0.6) for i in ax4.spines.itervalues()]
            
            for (y,yy) in zip(olat,mlat):

                if yy-y>0:

                    ax4.plot([y,y],[d2k(yy-y),0],'r-o',markevery=olat.size, ms=3, lw=.5, mew=.1, mec='k')

                else:

                    ax4.plot([y,y],[d2k(yy-y),0],'b-o',markevery=olat.size, ms=3, lw=.5, mew=.1, mec='k')

            ax4.set_xlabel('Latitude')
            ax4.set_ylabel('Dislocation (km)')
            ax4.set_ylim(-max([abs(_) for _ in ax4.get_ylim()]),max([abs(_) for _ in ax4.get_ylim()]))
            ax4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            xlim=ax4.get_xlim()
            ax4.plot(xlim,[0,0], linewidth=1, color='k')
            ax4.set_xlim(xlim)
            ax4.locator_params(axis='x',nbins=5)
            ax4.locator_params(axis='y',nbins=5)
            ax4.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
            
            ax5 = fig.add_subplot(gs[2:,2])
            [i.set_linewidth(0.6) for i in ax5.spines.itervalues()]

            for (x,y) in zip(odep,mdep):

                if y-x>0:

                    ax5.plot([x,x],[y-x,0],'r-o',markevery=odep.size, ms=3, lw=.5, mew=.1, mec='k')
                else:

                    ax5.plot([x,x],[y-x,0],'b-o',markevery=odep.size, ms=3, lw=.5, mew=.1, mec='k')

            ax5.set_xlabel('Depth (km)')
            ax5.set_ylabel('Dislocation (km)')
            ax5.set_ylim(-max([abs(_) for _ in ax5.get_ylim()]),max([abs(_) for _ in ax5.get_ylim()]))
            ax5.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            xlim=(0,self.dep_max)
            ax5.plot(xlim,[0,0], linewidth=1, color='k')
            ax5.set_xlim(xlim)
            ax5.locator_params(axis='x',nbins=5)
            ax5.locator_params(axis='y',nbins=5)
            ax5.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
            
            ax6 = fig.add_subplot(gs[2:,3])
            [i.set_linewidth(0.6) for i in ax6.spines.itervalues()]
            
            bns = arange(0, self.max_rms + self.max_rms/20.0, self.max_rms/20.0)
            ax6.hist([orms,mrms],color=['r','b'],label=['Initial','Relocated'], histtype='bar',bins=bns, rwidth=0.8, edgecolor='gray')
            ax6.set_xlabel('RMS (s)')
            ax6.set_ylabel('Event (#)')
            ax6.set_xlim(0,self.max_rms)
            ax6.locator_params(axis='x',nbins=4)
            ax6.locator_params(axis='y',nbins=5)
            ax6.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
            ax6.legend(loc=1, fontsize=5)

            print '+++ Minimum model id:',self.ind+1
            
            plt.tight_layout()
            plt.savefig(path.join('figs',self.project_nm,'velest_final_'+self.project_nm+'.tiff'))
            plt.close()
         
        # if no synthetic velocity model has been defined (single mode) * NOT IMPLEMENTED YET!

        else:

            pass

    #*-*-*-*-*- Plot data figures
    #
    # Plot some statistical figures
    # about data.
    
    def plot_data(self):

        inp_syntvel    = path.join('par','syntvel.dat')
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

        if not path.exists('figs'):

            mkdir('figs')

        if path.exists(path.join('figs',self.project_nm)):

            rmtree(path.join('figs',self.project_nm))
            mkdir(path.join('figs',self.project_nm))
    
        else:

            mkdir(path.join('figs',self.project_nm))

        cnv_dic         = ReadCNV(self.input_dic['EQ_datafile'])
        vel_dic,sta_dic = self.read_VelSta(inp_vel=path.join('velinp',self.velmod_name),inp_sta=self.input_dic['Stationfile'])

        self.sta_list = {}

        with open(path.join('tools','STATION0.HYP')) as f:

            for l in f:

                if (l[13:14] == 'N' or l[13:14] == 'S')  and (l[22:23] == 'E' or l[22:23] == 'W'):

                    n   = l[2:6].strip()
                    lat = ll.Latitude(degree=float(l[6:8]), minute=float(l[8:13])).decimal_degree
                    lon = ll.Longitude(degree=float(l[14:17]), minute=float(l[17:22])).decimal_degree
                    
                    self.sta_list[n] = {'lon':lon,'lat':lat}
                    
        print "+++ Plotting data statistics ..."

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

        self.orig_dep = evdep

        for _,__ in zip(range(len(vel_dic['P']['p_vel'])),range(len(vel_dic['P']['p_dep']))):

            vel.append(vel_dic['P']['p_vel'][_])
            vel.append(vel_dic['P']['p_vel'][__])
            dep.append(vel_dic['P']['p_dep'][__])
            dep.append(vel_dic['P']['p_dep'][__])
            
        dep.pop(0)
        dep.append(self.dep_max)

        # plot-No 1, velocity model

        init_plotting_isi(18,10)
        plt.rcParams['axes.labelsize'] = 7.5

        ax1 = plt.subplot(2,2,1)
        ax1.plot(vel,dep,color='black',linewidth=1.5,label='Refrence')
        xlim = ax1.get_xlim()
        ylim = ax1.get_ylim()
        xlim = (self.vel_min,self.vel_max)
        ylim = [min(dep),max(dep)+0.1*diff(ylim)[0]]
        ax1.set_xlim(xlim)
        ax1.set_ylim(ylim)
        ax1.set_ylim(ax1.get_ylim()[::-1])
        ax1.set_xlabel('Velocity (km/s)')
        ax1.set_ylabel('Depth (km)')
        ax1.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
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
        
        ax2.text(51.3890, 35.6892, '$Tehran$', clip_on=True)
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

        self.m.drawparallels(parallels,labels=[True,False,False,False],fmt='%4.1f',linewidth=0, fontsize=6)
        self.m.drawmeridians(meridians,labels=[False,False,False,True],fmt='%4.1f',linewidth=0, fontsize=6)

        #____ Scale & North sing

        self.scale_len = self.sc_len # km
        self.sc_lon_st = min(ax2.get_xlim())+0.04*diff(ax2.get_xlim())[0]
        self.sc_lon_ed = self.sc_lon_st + k2d(self.scale_len)
        self.sc_lat_st = min(ax2.get_ylim())+0.04*diff(ax2.get_ylim())[0]
        self.sc_lat_ed = self.sc_lat_st
        self.txt_x     = mean([self.sc_lon_st, self.sc_lon_ed])
        self.txt_y     = self.sc_lat_ed + abs(diff([self.sc_lat_ed,min(ax2.get_ylim())]))

        ax2.plot([self.sc_lon_st,self.sc_lon_ed],
                [self.sc_lat_st,self.sc_lat_ed],
                color='k',linewidth=1.5)
        ax2.plot(self.sc_lon_st,self.sc_lat_ed,
                color='k',linewidth=1.5, marker='>',ms=6)
        ax2.plot(self.sc_lon_ed,self.sc_lat_ed,
                color='k',linewidth=1.5, marker='<',ms=6)
        ax2.text(self.txt_x,self.txt_y,str(self.scale_len)+' km',
                horizontalalignment='center')

        ax2.annotate('N', xy=(.031, .07), xycoords='axes fraction', fontsize=10,  ha='center', va='bottom')
        ax2.annotate('^', xy=(.032, .09), xycoords='axes fraction', fontsize=15,  ha='center', va='bottom')

        dy = ax2.get_ylim()
        dy = dy[1] - dy[0]
        dy = dy
        fs     = 0.8*plt.rcParams['font.size']
        c      = evrms
        
        im  = ax2.scatter(evlon,evlat,c=c,marker='.',s=50,linewidth=None,edgecolor='black',cmap=plt.cm.jet, vmin=0.0, vmax=self.max_rms)
        ax2.plot(stlon,stlat, marker='^', ms=5, mec='k', mew=1.5, mfc='gray', ls='')
        ax2.plot(stlon[refid],stlat[refid], marker='^', ms=5, mec='gray', mew=1.5, mfc='k', ls='',label='Reference')
       
        for x,y,t in zip (stlon,stlat,stnm):

            if self.plt_stn == 'False': break
            
            if self.lon_min<=x<=self.lon_max and self.lat_min<=y<=self.lat_max:
                
                ax2.text(x,y-dy*.05,t,fontsize=fs, ha='center',bbox=dict(facecolor='w', edgecolor='k', boxstyle='round,pad=.2'))
                
        ax2.locator_params(axis='x',nbins=5)
        ax2.locator_params(axis='y',nbins=5)
        ax2.legend(loc=1)
        divider      = make_axes_locatable(ax2)
        cax          = divider.append_axes("right", size="4%", pad=0.05)
        cb           = plt.colorbar(im, ax=ax2, cax=cax, label='RMS (s)')
        tick_locator = ticker.MaxNLocator(nbins=6)
        cb.locator   = tick_locator
        cb.update_ticks()
        ax2.set_xlim(self.lon_min,self.lon_max)
        ax2.set_ylim(self.lat_min,self.lat_max)
        
        # plot-No 3, depth-number_of_station distribution

        ax3 = plt.subplot(2,2,3)
        c   = evrms
        im  = ax3.scatter(evdep,nump,c=c,marker='.',s=100,edgecolor='black',cmap=plt.cm.jet, vmin=0.0, vmax=self.max_rms)
        ax3.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
        ax3.set_xlabel('Depth (km)')
        ax3.set_ylabel('Total Number of Phase (P,S)')
        ax3.set_xlim(0,max(evdep))
        ax3.set_ylim(0,max(ax3.get_ylim()))
        ax3.locator_params(axis='x',nbins=5)
        ax3.locator_params(axis='y',nbins=5)
        divider      = make_axes_locatable(ax3)
        cax          = divider.append_axes("right", size="4%", pad=0.05)
        cb           = plt.colorbar(im, ax=ax3, cax=cax, label='RMS (s)')
        tick_locator = ticker.MaxNLocator(nbins=6)
        cb.locator   = tick_locator
        cb.update_ticks()

        # plot-No 4, gap-number_of_station distribution

        ax4 = plt.subplot(2,2,4)
        c   = evrms
        im  = ax4.scatter(evgap,nump,c=c,marker='.',s=100,edgecolor='black',cmap=plt.cm.jet, vmin=0.0, vmax=self.max_rms)
        ax4.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
        ax4.set_xlabel('Azimuthal Gap ($\circ$)')
        ax4.set_ylabel('Total number of Phase (P,S)')
        ax4.set_xlim(min(evgap),max(evgap))
        ax4.set_ylim(0,max(ax4.get_ylim()))
        ax4.locator_params(axis='x',nbins=5)
        ax4.locator_params(axis='y',nbins=5)
        divider      = make_axes_locatable(ax4)
        cax          = divider.append_axes("right", size="4%", pad=0.05)
        cb           = plt.colorbar(im, ax=ax4, cax=cax, label='RMS (s)')
        tick_locator = ticker.MaxNLocator(nbins=6)
        cb.locator   = tick_locator
        cb.update_ticks()
        
        plt.tight_layout()
        plt.savefig(path.join('figs',self.project_nm,'data_'+self.project_nm+'.tiff'))
        plt.close()


    ############################################## 
    #*-*-*-*-*- Make final reports
    #
    # Make final reports
    # them.
    
    def mk_report(self):

        print "+++ Preparing reports ..."

        self.repfi    = open(path.join('figs',self.project_nm,'report.dat'),'w')
        self.modout   = open(path.join('figs',self.project_nm,'model.mod'),'w')

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
            
            self.repfi.write(' %4.2f     %4.1f\n'%(v,d))

        self.repfi.write('\n# 2- Mean velocity model:\n\n')

        self.repfi.write('Vel[km/s]   Dep[km]\n')

        self.mean_mod[1][0] = 0.0 # set first depth to zero
        
        for v,d in zip(self.mean_mod[0][::2],self.mean_mod[1][::2]):
            
            sel
