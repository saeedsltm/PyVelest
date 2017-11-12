import matplotlib.pyplot as plt

def init_plotting():

    plt.rcParams['figure.figsize']       = (16, 9)
    plt.rcParams['figure.dpi']           = 75
    plt.rcParams['font.size']            = 16
    plt.rcParams['font.family']          = 'Times New Roman'
    plt.rcParams['axes.labelsize']       = 18
    plt.rcParams['axes.titlesize']       = 20
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
    plt.rcParams['savefig.dpi']          = 300
    plt.rcParams['xtick.minor.visible']  = 'False'
    plt.rcParams['ytick.minor.visible']  = 'False'
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')
    plt.locator_params(nticks=4)

def cm2inch(value):
    return value/2.54

def init_plotting_isi(*args):

    if len(args) == 2:

        xsize = cm2inch(float(args[0]))
        ysize = cm2inch(float(args[1]))

    else:

        xsize = cm2inch(9.5)
        ysize = cm2inch(11.)

    plt.rcParams['figure.figsize']       = (xsize,ysize)
    plt.rcParams['figure.dpi']           = 300
    plt.rcParams['font.size']            = 8
    plt.rcParams['font.family']          = 'Times New Roman'
    plt.rcParams['axes.labelsize']       = 7
    plt.rcParams['axes.titlesize']       = 8
    plt.rcParams['legend.fontsize']      = 7
    plt.rcParams['xtick.labelsize']      = 7
    plt.rcParams['ytick.labelsize']      = 7
    plt.rcParams['xtick.major.size']     = 2
    plt.rcParams['xtick.minor.size']     = 1.5
    plt.rcParams['xtick.major.width']    = 0.50
    plt.rcParams['xtick.minor.width']    = 0.25
    plt.rcParams['ytick.major.size']     = 2
    plt.rcParams['ytick.minor.size']     = 1.5
    plt.rcParams['ytick.major.width']    = 0.50
    plt.rcParams['ytick.minor.width']    = 0.25
    plt.rcParams['legend.frameon']       = True
    plt.rcParams['legend.shadow']        = False
    plt.rcParams['legend.loc']           = 'lower left'
    plt.rcParams['legend.numpoints']     = 1
    plt.rcParams['legend.scatterpoints'] = 1
    plt.rcParams['axes.linewidth']       = .5
    plt.rcParams['savefig.dpi']          = 300
    plt.rcParams['xtick.minor.visible']  = 'False'
    plt.rcParams['ytick.minor.visible']  = 'False'
    plt.rcParams['patch.linewidth']      = .5
