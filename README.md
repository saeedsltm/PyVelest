Short manual for "How To Run" PyVelest/Pycheck and other modules.

Requirements:

  *NOTE!*

*** It assumed that you are familiar with the VELEST program already, knowing how it works, how to run it and how to 
        interpret its results. Otherwise this code is not what you're looking for. 

*** Before using the code, please be sure that the following programs are already installed on your Linux system (-
        tested on Ubuntu 16.04).

Non-Pythonic program:

- SEISAN, (http://seisan.info/)

Pythonic programs:

- Numpy, (http://www.numpy.org/)
- Scipy, (https://www.scipy.org/)
- Matplotlib, (https://matplotlib.org/)
- Statsmodels, (http://www.statsmodels.org/stable/index.html)
- LatLon (https://pypi.python.org/pypi/LatLon)

Input files:

    There is only needed to prepare two input files. One is phase arrival time with NORDIC format (e.g,. 
        collect.out) and the other one should include station information like STATION0.HYP file used by SEISAN (look at
        SEISAN manual).


Processing:

    After two input files have been prepared, put them into "tools" directory, then open a terminal in "tools" and do the 
       followings:

    Note: root directory = the directory that you unzipped the package.

    1- Select best event: Use modules "PySelect2d.py" for 2D selection or "PySelect3d.py" for 3D selection. These modules
       will read the input file "select2D.par" or "select3D.par" depends  what you've selected. Please set the parameters
       inside the input parameter file according to your region and based on the criteria you want to select best events.
       The input parameters file are self-explanatory. Then (here for 2D case) type "python PySelect2D.py"  and  give the 
       input NORDIC file name. Check if the output file "select.out" is created.

![select3d](https://user-images.githubusercontent.com/20725501/32775311-5eac46c2-c92f-11e7-9b5f-d0a2e283129a.png)

    2- Check which station could be removed. If number of recorded phases in any station is less than let say 10, and it
       shows high RMS, i recommend to remove that station from your station list "STATION0.HYP". To check if there is any
       such a station, type "python PyNordicStat.py" and give the NORDIC input file name. It generate a report file called
       "statistic_report.out" and a figure called "statistical_plot.tiff". Check these files and decide if it is necessary
       to remove any or not.

    3- Plot Vp/Vs. If you are not sure about Vp/Vs ratio, you can calculate this value using "PyVpVs.py". To do that, type
       "python PyVpVs.py" in terminal and give the NORDIC input file name. It generate a figure called "VpVs_x.tiff" which
        "x" is the Vp/Vs ratio. This module calculats the ratio aftre removing outliers from your dataset.

![vpvs_1 73](https://user-images.githubusercontent.com/20725501/32775531-ef81b150-c92f-11e7-8c57-eb3e0ed4644f.png)

    4- Generate VELEST input file automatically. Using the module "PyVelestInp.py" it creates all input files required by 
       VELEST program. Type "python PyVelestInp.py" and give the NORDIC input file name. it generates the "velinp" folder
       in root, and phase data in CNV format "data.cnv", velocity model "initial.mod" and station file "station.sta". the
       reference station will be selected automatically by module in terms of a) number of recorded phase, b) RMS residual
       , mean distance from the center of region study.

    5- You may now close the terminal and go back to root directory.

    6- Now you must set the required parameters by VESTEST for inversion, making synthetic models and finally for plotting
       results. Go to "par" directory where you find three files inside. Edit the file "par.dat" and set only the required
       parameters for running VELEST. Then you should edit the file "syntvel.dat" for make synthetic models. Remember to
       choose different name for 'velmod_name' parameter, each time you want to do a new run. Finally set the parameters on
       "topo.dat" file which will be used for plotting figures.

    7- Open a terminal in root directory then type "python PyVelest.py". Given the name you already set in  "syntvel.dat"
       file for 'velmod_name' parameter and hit the 'Enter' to do the first run. For each new run (after you modified all
       parameter files or initial dataset/model) give a new name and run. At the end of each run, a new directory will be 
       created inside "figs" folder. The following is a short discription of each file:

       a) bad_select.out : A NORDIC phase file containing all bad events.
       b) fin_hyp.cnv : Final phase file in CNV format.
       c) model.mod : Best model obtained after inversion.
       d) report.dat : A summary report file.
       e) select.out : Final phase file in NORDIC format.
       f) sta_cor.out : station corrections obtained after inversion with best fitted model.
       h) data_initial.tiff : Figure indicating initial data status.
       i) final_initial.tiff : Figure indicating final results.
       j) path_initial.tiff : Figure indicating path coverage both in horizontal and depth.
       k) synt_models_initial.tiff : Figure indicating generated synthetic models.
       l) topo_diff_initial.tiff : Figure indicate epicenter dislocation between initial and final locations using vectors.
       m) topo_fin.tiff : Figure indicating final epicenters with error ellipse.
       n) topo_ini.tiff : Figure indicating initial epicenters with error ellipse.
       o) velest_analytics_initial.tiff : Figure indicating VELEST analyzing in each iteration step.
       p) velest_final_initial.tiff : Figure indicating final results of VELEST outputs.

![data_initial](https://user-images.githubusercontent.com/20725501/32776365-9c1862cc-c932-11e7-9c91-6ee7a54597c1.png)
![path_initial](https://user-images.githubusercontent.com/20725501/32776414-c6656110-c932-11e7-948e-c9ac23852f93.png)
![synt_models_initial](https://user-images.githubusercontent.com/20725501/32776450-e123ddc4-c932-11e7-946c-9c8a0f7a973d.png)
![velest_analytics_initial](https://user-images.githubusercontent.com/20725501/32776491-045c8bb0-c933-11e7-976c-da98c55be26f.png)
![final_initial](https://user-images.githubusercontent.com/20725501/32776530-24740e28-c933-11e7-90cf-fae4da8bdcc9.png)

    8- After you finish with the "PyVelest.py" code, it's time to check the stability of the model and make some comparison
       between calculated models. In some cases you may also to know about the best  damping  values. To do these, there is
       "PyCheck.py" module in root directory. Running this code is similar to "PyVelest.py"  so  you will be asked in each
       step. The output figures of running this code, will be saved inside "Check-Test" folder and its subdirectories.

![velocity_damping_test](https://user-images.githubusercontent.com/20725501/32776663-94e9c742-c933-11e7-91bb-b2e715859003.png)
![stabilitytest](https://user-images.githubusercontent.com/20725501/32776683-a40d6332-c933-11e7-8ab5-6de31bce9fa7.png)

The Code; Files and Folders:

    Inside root directory, you will find 3 folders and 3 files (including me). The following is  a  short description about
    each of them:
    
    Folders:
           - gdb : This folder includes all files and folders for topography maps currently support for Iran rigion. So in
                   this version there is no support to plot topography maps outside of Iran. You  may  Not consider it but,
                   you should not remove or modify its component. Otherwise the code may not work.    
           - par : This folder includes all parameter files required for running VELEST, making synthetic models and finally
                   plotting figure.
           - tools : This folder and its components/sub-folders are combination of codes, input data and parameter files.

    Files:
         - tools:
                - utility : This folder contains only some python codes to read or write data formats and plotting defaults.
                - addnoise.par : input parameter file used for generating noisy dataset. This will be used when you want to
                                 do stability test. 
                - PyAddNoise.py : Module for generating noisy dataset. Before running this code  you  must modify the input 
                                  parameter file "addnoise.par".
                - PyNordicStat.py : Make a report and figure describing initial status of used dataset.
                - PySelect2D.py : This module is used for make a robust selection of events in 2D.
                - PySelect3D.py : This module is used for make a robust selection of events in 3D.
                - PyVelestInp.py : This module is used for preparing VELEST input files.
                - PyVpVs.py : This module is used for calculating the Vp/Vs ratio.
                - select2D.par : Input parameter file required by PySelect2D.py.
                - select3D.par : Input parameter file required by PySelect3D.py.

GOOD-LUCK!
