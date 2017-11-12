import pylab as plt
from glob import glob
import os, sys
from shutil import copy
from numpy.linalg import norm
from scipy.stats import gaussian_kde
from numpy import arange, array, pi, sqrt, degrees, arctan2, linalg, loadtxt, mgrid, reshape, vstack, hstack, linspace
from mpl_toolkits.axes_grid1 import make_axes_locatable
from copy import deepcopy as dp
sys.path.append(os.path.join('tools','utility'))
from initial_mpl import init_plotting_isi
from matplotlib import ticker, gridspec

#+++++++++++++++++++ PLOTTING VELOCIT MODEL SECTION

def k2d(kilometer, radius=6371):

    return kilometer / (2.0 * radius * pi / 360.0)


def d2k(degrees, radius=6371):

    return degrees * (2.0 * radius * pi / 360.0)

#___________________ EXTRACT INITIAL RESULT

def f_ini(mod_file):

    d = []
    v = []
    c = 0 

    with open(mod_file) as f:

        start_flag = False

        for l in f:

            if 'P-VELOCITY MODEL' in l:

                start_flag = True

            if start_flag and len(l.split()) == 1:

                start_flag = False

            if start_flag:

                v.append(float(l.split()[0]))
                d.append(float(l.split()[1]))

            c+=1
            
    return v,d    

#___________________ EXTRACT FINAL RESULT

def f_fin(rep_file):
    
    d = []
    v = []

    with open(rep_file) as f:

        end_flag = False
        start_flag = False

        for l in f:
            

            if 'Vel[km/s]   Dep[km]' in l:

                start_flag = True
                
            if start_flag and not l.split():
                
                end_flag = True
                start_flag = False

            if start_flag and end_flag and 'Vel[km/s]   Dep[km]' not in l:

                v.append(float(l.split()[0]))
                d.append(float(l.split()[1]))

    return v,d

#___________________ PLOT VELOCITY MODELS

def plot_vel():

    init_plotting_isi(8,8)
    plt.rcParams['axes.labelsize'] = 7

    max_dep = input('\n+++ Maximum depth to plot:\n')
    vel_rng = input('\n+++ Minimum, Maximum velocity range to plot:\n')
    
    velmods  = sorted(glob(os.path.join('velinp','*.mod')))
    projects = sorted(glob(os.path.join('figs','*')))
    color    = ['b','g','r','c','m','y','k',]

    for velmod, project, c in zip(velmods, projects, color):

        rep_file_ini = os.path.join(velmod)
        rep_file_fin = os.path.join(project,'report.dat')

        v_ini,d_ini  = f_ini(rep_file_ini)
        title_ini    = rep_file_ini.split(os.sep)[-1].split('.')[0]
        d_ini        = [-1*_ for _ in d_ini]
        
        v_fin,d_fin  = f_fin(rep_file_fin)
        title_fin    = rep_file_fin.split(os.sep)[1]
        d_fin        = [-1*_ for _ in d_fin]

        x_ini  = list(hstack([[m,n] for m,n in zip(v_ini,v_ini)]))
        y_ini  = list(hstack([[m,n] for m,n in zip(d_ini,d_ini)]))

        x_fin  = list(hstack([[m,n] for m,n in zip(v_fin,v_fin)]))
        y_fin  = list(hstack([[m,n] for m,n in zip(d_fin,d_fin)]))

        y_ini.pop(0)
        y_ini.append(-1*max_dep)
        
        y_fin.pop(0)
        y_fin.append(-1*max_dep)        
            
        plt.plot(x_fin, -array(y_fin), linewidth=1,label=title_fin,color=c, zorder=10,alpha=.9)
        plt.plot(x_ini, -array(y_ini), linewidth=1,label=title_ini,color=c, linestyle='--',zorder=10,alpha=.9)
        
        plt.ylim(ymin=max_dep,ymax=-0.1)
        plt.xlim(vel_rng[0],vel_rng[1])
        plt.xlabel("Velocity (km/s)")
        plt.ylabel("Depth (km)")
        plt.legend(loc=1, fontsize=6)
        plt.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
        [i.set_linewidth(0.6) for i in plt.gca().spines.itervalues()]

    plt.tight_layout()    
    plt.savefig(os.path.join('ModelComparsion.tiff'),dpi=300)
    plt.close()

    print '\n+++ Output file "ModelComparsion.tiff" was created.\n'


#+++++++++++++++++++ DAMPPING TESTS SECTION


#___________________ EXTRACT FINAL RESULTS

def read_velestout(inp='velest.out'):

    # READ BEST VELEST INVERSION RESULT 

    model_var = []
    data_var  = []
    
    with open(inp) as f:

        flag = False
        hint = ' Velocity model   1'

        for l in f:

            if hint in l:

                flag = True
                tmp  = []
                
            if flag and not l.strip():

                model_var.append(norm(tmp)**2)
                flag = False

            if flag and hint not in l:

                tmp.append(float(l.split()[1]))

    with open(inp) as f:

        hint = 'DATVAR='
        
        for l in f:

            if hint in l:

                    data_var.append(float(l.split()[1]))

    data_var.pop(0) 
    return data_var, model_var
               
        
#___________________ DAMPING TEST

def damping_test():
    
    project_nm = raw_input('\n+++ Enter project name:\n\n')
    velmod_nm  = raw_input('\n+++ Enter best inversion id number (e.g. minimum model id):\n\n')
    which_dmp  = input('\n+++ Which damping parameter do you want to test?\n\n1- othet (Origin Time Damping Factor)\n2- xythet (Epicenter Damping Factor)\n3- zthet (Depth Damping Factor)\n4- vthet (Velocity Damping Factor)\n5- stathet (Station correction Damping Factor)\n\n')
    vel_dmp    = input('\n+++ selected damping range (Min, Max, Num of damping point):\n\n')

    if not os.path.exists(os.path.join('Check-Test',project_nm)):

        os.makedirs(os.path.join('Check-Test',project_nm))
        
    damp_dic = {1:'Origin Time Damping Test',
                2:'Epicenter Damping Test',
                3:'Depth Damping Test',
                4:'Velocity Damping Test',
                5:'Station Correction Damping Test',}

    copy(os.path.join('velout',project_nm,'velest%s.out'%(velmod_nm)), 'velest.cmn')

    print '\n+++ Warning! selected dampping parameter for the following file will be changed!'
    print '\n%s\n'%(os.path.join('velout',project_nm,'velest%s.out'%(velmod_nm)))

    tmp = open('tmp.dat', 'w')
    tmp.write('damping test')

    with open('velest.cmn') as f:

        c = 0

        for l in f:

            if 31 <= c <= 129:

                tmp.write(l)

            c+=1

    tmp.close()

    model_var = []
    data_var  = []

    vel_dmp = linspace(vel_dmp[0], vel_dmp[1], vel_dmp[2], dtype=int)

    flag = False

    print ''

    for vd in vel_dmp:

        print '+++ Test on Damping = %4d'%(vd)

        output = open('velest.cmn','w')

        with open('tmp.dat') as f:

            hint = '***   othet   xythet    zthet    vthet   stathet'

            for l in f:

                if hint in l:

                    flag = True

                if flag and hint not in l:

                    l    = l.split()
                    l[which_dmp-1] = '%.3f'%(vd)
                    l    = '     '+'    '.join(l)+'\n'
                    flag = False

                output.write(l)
                
        output.close()

        os.system('velest > /dev/null')

        x, y = read_velestout(inp=os.path.join('velout',project_nm,'velest%s.out'%(velmod_nm)))

        data_var.append(x)
        model_var.append(y)

    data_var  = array(data_var)
    model_var = array(model_var)

    init_plotting_isi(8,8)
    plt.rcParams['axes.labelsize'] = 7
    
    ax = plt.subplot(111)
    [i.set_linewidth(0.6) for i in ax.spines.itervalues()]
    
    ax.plot(model_var[:,0],data_var[:,0],marker='o', ms=3, lw=1, color='r')
    for x,y,z in zip(model_var[:,0], data_var[:,0], vel_dmp): ax.text(x,y,'%d'%(z))
    ax.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
    ax.set_title('%s'%(damp_dic[which_dmp]))
    ax.set_xlabel('Model Variance (km/s)$^2$')
    ax.set_ylabel('Data Variance (s)$^2$')

    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.locator_params(axis='x',nbins=5)
    ax.locator_params(axis='y',nbins=5)
            
    name = damp_dic[which_dmp].replace(' ','_')

    plt.tight_layout()
    plt.savefig(os.path.join('Check-Test',project_nm,name+'.tiff'),dpi=300)
    plt.close()

    print '\n+++ Output file "%s" was created.\n'%(os.path.join('Check-Test',project_nm,name+'.tiff'))


#+++++++++++++++++++ STABILITY TEST ON FINAL RESULT


def stability_test():
    
    project_name = raw_input('\n+++ Enter project name (project that returned the best results):\n\n')
    model_mod    = raw_input('\n+++ Which velocity model should to be used?\n\n1- Mean Model.\n2- Minimum Model (default).\n\n')
    max_dep      = raw_input('\n+++ Maximum depth to plot velocity model (def=30.0km):\n\n')

    print '\n+++ Please set the following damping and inversion mode parameters:\n\n'

    othet   = raw_input('Origin Time (def=0.01)       = ')
    xythet  = raw_input('Epicenter (def=0.01)         = ')
    zthet   = raw_input('Depth (def=0.01)             = ')
    vthet   = raw_input('Velocity (def=1.0)           = ')
    stathet = raw_input('Station Correction (def=1.0) = ')
    invrat  = raw_input('Invert ratio (def=2)         = ')
    numitr  = raw_input('Number of iteration (def=7)  = ')

    dmp_def_val = {'othet':0.01, 'xythet':0.01, 'zthet':0.01, 'vthet':1.0, 'stathet':1.0, 'invrat':2, 'numitr':7, 'max_dep':30.0}

    for i in dmp_def_val.keys():

        dmp = eval(i)

        if dmp:

            dmp_def_val[i] = float(dmp)

    if not os.path.exists(os.path.join('Check-Test',project_name)):

        os.makedirs(os.path.join('Check-Test',project_name))

    # PREPARE REQUIRED FILES

    vel_mod_file = open(os.path.join('Check-Test',project_name,'model.mod'),'w')

    with open(os.path.join('figs',project_name,'model.mod')) as f:

        flag = False

        for l in f:

            if model_mod == '1':

                if 'Model: mean.mod' in l: flag = True
                if flag and not l.strip(): flag = False

                if flag:

                    vel_mod_file.write(l)

            else:

                if 'Model: min.mod' in l: flag = True
                if flag and not l.strip(): flag = False

                if flag:

                    vel_mod_file.write(l)            

    vel_mod_file.close()

    with open(os.path.join('figs',project_name,'report.dat')) as f:

        for l in f:

            if 'ID number' in l:

                ind = float(l.split()[-1])
                
    copy(os.path.join('velout',project_name,'velest%d.out'%(ind)), 'tmp.dat')

    velestcmn = open(os.path.join('Check-Test',project_name,'velest.cmn'), 'w')
    velestcmn.write('Velocity Stability Test')

    with open('tmp.dat') as f:
        
        h0 = '***   othet   xythet    zthet    vthet   stathet'
        h1 = '*** Modelfile:'
        h2 = '*** Stationfile:'
        h3 = '*** File with Earthquake data:'
        h4 = '*** Main print output file:'
        h5 = '*** File with final hypocenters in *.cnv format:'
        h6 = '*** File with new station corrections:'
        h7 = '*** delmin   ittmax   invertratio'

        f0         = False 
        f1, f2, f3 = False, False, False
        f4, f5, f6 = False, False, False
        f7         = False

        c = 0

        for l in f:

            if 31 <= c <= 121:

                if h0 in l:

                    f0 = True

                if f0 and h0 not in l:

                    velestcmn.write('     %.3f     %.3f     %.3f     %.3f     %.3f\n'%(dmp_def_val['othet'],
                                                                                       dmp_def_val['xythet'],
                                                                                       dmp_def_val['zthet'],
                                                                                       dmp_def_val['vthet'],
                                                                                       dmp_def_val['stathet']))
                    f0 = False
                    continue
                
                if h1 in l:

                    f1 = True

                if f1 and h1 not in l:

                    velestcmn.write('model.mod\n')
                    f1 = False
                    continue


                if h2 in l:

                    f2 = True

                if f2 and h2 not in l:

                    velestcmn.write(os.path.join('..'+os.sep+'..','figs',project_name,'sta_cor.out')+'\n')
                    f2 = False
                    continue


                if h3 in l:

                    f3 = True

                if f3 and h3 not in l:

                    velestcmn.write(os.path.join('..'+os.sep+'..','velinp','noisy.cnv')+'\n')
                    f3 = False
                    continue


                if h4 in l:

                    f4 = True

                if f4 and h4 not in l:

                    velestcmn.write('velest.out\n')
                    f4 = False
                    continue


                if h5 in l:

                    f5 = True

                if f5 and h5 not in l:

                    velestcmn.write('final_loc.cnv\n')
                    f5 = False
                    continue


                if h6 in l:

                    f6 = True

                if f6 and h6 not in l:

                    velestcmn.write('station_cor.sta\n')
                    f6 = False
                    continue


                if h7 in l:

                    f7 = True

                if f7 and h7 not in l:

                    velestcmn.write('    0.010       %d          %d\n'%(dmp_def_val['numitr'],
                                                                        dmp_def_val['invrat']))
                    f7 = False
                    continue
              
                velestcmn.write(l)

            c+=1

    velestcmn.close()

    # RUN VELEST TO PERFORM TEST

    here = os.getcwd()
    os.chdir(os.path.join('Check-Test',project_name))
    os.system('velest > /dev/null')
    os.chdir(here)

    # PLOT STABILITY RESULTS

    x_ini = []
    y_ini = []
    z_ini = []

    x_nsy = []
    y_nsy = []
    z_nsy = []

    x_fin = []
    y_fin = []
    z_fin = []

    # DEFINE RESOURCES

    d_ini = os.path.join('figs',project_name,'fin_hyp.cnv')
    d_nsy = os.path.join('velinp','noisy.cnv')
    d_fin = os.path.join('Check-Test',project_name,'final_loc.cnv')

    for x,y,z,inp in zip([x_ini,x_nsy,x_fin],[y_ini,y_nsy,y_fin],[z_ini,z_nsy,z_fin],[d_ini,d_nsy,d_fin]):
        
        with open(inp) as f:

            for l in f:

                if l[25:26] == 'N' and l[35:36] == 'E':

                    y.append(float(l[17:25]))
                    x.append(float(l[26:35]))
                    z.append(float(l[36:43]))

    vel_ini = [[],[]]

    with open(os.path.join('Check-Test',project_name,'model.mod')) as f:

        flag = False

        for l in f:

            if 'P-VELOCITY MODEL' in l: flag = True
            if flag and len(l.split()) == 1: break

            if flag:

                vel_ini[0].append(float(l.split()[0]))
                vel_ini[1].append(float(l.split()[1]))

    vi = dp(vel_ini[0])
    y  = [-1*d for d in vel_ini[1]]
    x  = list(hstack([[m,n] for m,n in zip(vel_ini[0],vel_ini[0])]))
    y  = list(hstack([[m,n] for m,n in zip(y,y)]))
    y.pop(0)
    y.append(-1*dmp_def_val['max_dep'])        

    vel_ini = [x,y]

    with open(os.path.join('Check-Test',project_name,'velest.out')) as f:

        flag = False

        for l in f:

            if ' Velocity model   1' in l:

                flag    = True
                vel_fin = [[],[]]
                continue

            if flag and not l.strip(): flag = False

            if flag and l.split()[0][0].isdigit():

                vel_fin[0].append(float(l.split()[0]))
                vel_fin[1].append(float(l.split()[2]))

    vf = dp(vel_fin[0])
    y  = [-1*d for d in vel_fin[1]]
    x  = list(hstack([[m,n] for m,n in zip(vel_fin[0],vel_fin[0])]))
    y  = list(hstack([[m,n] for m,n in zip(y,y)]))
    y.pop(0)
    y.append(-1*dmp_def_val['max_dep'])        

    vel_fin = [x,y]

    init_plotting_isi(17,9)
    plt.rcParams['axes.labelsize'] = 7
    
    ax1 = plt.subplot(221)

    x_diff =  d2k(array(x_nsy) - array(x_ini))
    xma, xmi = max(x_ini), min(x_ini)
    yma, ymi = max(x_diff), min(x_diff)
    ax1.plot(x_ini, x_diff, color='r', marker='x', ms=3, linestyle='', zorder=102)
    x_diff =  d2k(array(x_fin) - array(x_ini))
    ax1.plot(x_ini, x_diff, color='b', marker='x', ms=3, linestyle='', zorder=102)
    ax1.set_xlim(xmi-(xma-xmi)*.05, xma+(xma-xmi)*.05)
    ax1.set_ylim(ymi-(yma-ymi)*.05, yma+(yma-ymi)*.05)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Dislocation (km)')
    ax1.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
    ax1.locator_params(axis='x',nbins=6)
    ax1.locator_params(axis='y',nbins=6)

    # Plot KDE

    xmin, xmax = ax1.get_xlim()
    ymin, ymax = ax1.get_ylim()  
    X, Y       = mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions  = vstack([X.ravel(), Y.ravel()])
    values     = vstack([x_ini, x_diff])
    kernel     = gaussian_kde(values)
    Z          = reshape(kernel(positions).T, X.shape)
    im         = ax1.contourf(X, Y, Z, cmap=plt.cm.gist_earth_r, alpha=.9, zorder=101)
    divider    = make_axes_locatable(ax1)
    cax        = divider.append_axes("right", size="4%", pad=0.05)
    cb         = plt.colorbar(im, ax=ax1, cax=cax)
    tk_locator = ticker.MaxNLocator(nbins=6)
    cb.locator = tk_locator
    cb.update_ticks()
    cb.outline.set_linewidth(.75)
    cb.set_label(label='PDF', size=6)
    cb.ax.tick_params(labelsize=6) 
          
    ax2 = plt.subplot(222)

    y_diff =  d2k(array(y_nsy) - array(y_ini))
    xma, xmi = max(y_ini), min(y_ini)
    yma, ymi = max(y_diff), min(y_diff)
    ax2.plot(y_ini, y_diff, color='r', marker='x', ms=3, linestyle='', zorder=102)
    y_diff =  d2k(array(y_fin) - array(y_ini))
    ax2.plot(y_ini, y_diff, color='b', marker='x', ms=3, linestyle='', zorder=102)
    ax2.set_xlim(xmi-(xma-xmi)*.05, xma+(xma-xmi)*.05)
    ax2.set_ylim(ymi-(yma-ymi)*.05, yma+(yma-ymi)*.05)
    ax2.set_xlabel('Latitude')
    ax2.set_ylabel('Dislocation (km)')
    ax2.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
    ax2.locator_params(axis='x',nbins=6)
    ax2.locator_params(axis='y',nbins=6)

    # Plot KDE

    xmin, xmax = ax2.get_xlim()
    ymin, ymax = ax2.get_ylim()
    X, Y       = mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions  = vstack([X.ravel(), Y.ravel()])
    values     = vstack([y_ini, y_diff])
    kernel     = gaussian_kde(values)
    Z          = reshape(kernel(positions).T, X.shape)
    im         = ax2.contourf(X, Y, Z, cmap=plt.cm.gist_earth_r, alpha=.9, zorder=101)
    divider    = make_axes_locatable(ax2)
    cax        = divider.append_axes("right", size="4%", pad=0.05)
    cb         = plt.colorbar(im, ax=ax2, cax=cax)
    tk_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tk_locator
    cb.update_ticks()
    cb.outline.set_linewidth(.75)
    cb.set_label(label='PDF', size=6)
    cb.ax.tick_params(labelsize=6)
    
    ax3 = plt.subplot(223)

    z_diff =  array(z_nsy) - array(z_ini)
    xma, xmi = max(z_ini), min(z_ini)
    yma, ymi = max(z_diff), min(z_diff)
    ax3.plot(z_ini, z_diff, color='r', marker='x', ms=3, linestyle='', zorder=102)
    z_diff =  array(z_fin) - array(z_ini)
    ax3.plot(z_ini, z_diff, color='b', marker='x', ms=3, linestyle='', zorder=102)
    ax3.set_xlim(xmi-(xma-xmi)*.05, xma+(xma-xmi)*.05)
    ax3.set_ylim(ymi-(yma-ymi)*.05, yma+(yma-ymi)*.05)
    ax3.set_xlabel('Depth (km)')
    ax3.set_ylabel('Dislocation (km)')
    ax3.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
    ax3.locator_params(axis='x',nbins=6)
    ax3.locator_params(axis='y',nbins=6)

    # Plot KDE

    xmin, xmax = ax3.get_xlim()
    ymin, ymax = ax3.get_ylim()
    X, Y       = mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions  = vstack([X.ravel(), Y.ravel()])
    values     = vstack([z_ini, z_diff])
    kernel     = gaussian_kde(values)
    Z          = reshape(kernel(positions).T, X.shape)
    im         = ax3.contourf(X, Y, Z, cmap=plt.cm.gist_earth_r, alpha=.9, zorder=101)
    divider    = make_axes_locatable(ax3)
    cax        = divider.append_axes("right", size="4%", pad=0.05)
    cb         = plt.colorbar(im, ax=ax3, cax=cax)
    tk_locator = ticker.MaxNLocator(nbins=6)
    cb.locator = tk_locator
    cb.update_ticks()
    cb.outline.set_linewidth(.75)
    cb.set_label(label='PDF', size=6)
    cb.ax.tick_params(labelsize=6)
    
    ax4 = plt.subplot(224)
    [i.set_linewidth(0.6) for i in ax4.spines.itervalues()]
    
    ax4.plot(vel_ini[0], -array(vel_ini[1]), 'r-', linewidth=1.5, label='Initial')
    ax4.plot(vel_fin[0], -array(vel_fin[1]), 'b-', linewidth=1.5, label='Inverted')
    model_rms_fit = norm(array(vf)-array(vi))
    ax4.set_xlabel('Velocity (km/s)')
    ax4.set_ylabel('Depth (km)')
    ax4.invert_yaxis()
    ax4.text(0.82, 0.9,'||model||=%.2f (km/s)'%(model_rms_fit), ha='center', va='center', transform=ax4.transAxes,
             bbox=dict(facecolor='w', alpha=0.5, pad=2), fontsize=6)
    ax4.grid(True, linestyle='--', linewidth=.5, color='k', alpha=.3)
    ax4.locator_params(axis='x',nbins=6)
    ax4.locator_params(axis='y',nbins=5)
    ax4.set_ylim(-min(vel_ini[1]),0)
    plt.legend(loc=3, fontsize=6)


    plt.tight_layout()
    plt.savefig(os.path.join('Check-Test',project_name,'StabilityTest.tiff'),dpi=300)
    plt.close()

    print '\n+++ Result was saved in "Check-Test%s" directory.\n'%(os.sep+project_name)

#___________________ RUN

ans = input('\n+++ Select one of the following options:\n\n1- Plot initial VS final velocity models.\n2- Perform dampping test.\n3- Perform stability test.\n\n')
    
if ans == 1:

    plot_vel()

if ans == 2:

    damping_test()

    for i in ['tmp.dat', 'velest.cmn']:

        if os.path.exists(i): os.remove(i)

if ans == 3:

    stability_test()

    for i in ['tmp.dat']:

        if os.path.exists(i): os.remove(i)
