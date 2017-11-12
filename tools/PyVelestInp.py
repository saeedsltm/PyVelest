import os
import sys
sys.path.append(os.path.join('utility'))
from shutil import rmtree, move
from PyNordicRW import Read_Nordic
from PyNordic2CNV import Nordic2CNV
from PyGetRefSta import PRD

nordic_inp = raw_input('\nEnter NORDIC file name:\n\n')

print ''

velmenu_inp = open('velmenu.inp','w')

velmenu_inp.write(nordic_inp+'\n')
velmenu_inp.write('1\n')
velmenu_inp.write('I\n')
velmenu_inp.write('3\n')
velmenu_inp.write('5\n')
velmenu_inp.write('A\n')
velmenu_inp.write('n\n')
velmenu_inp.write('Q\n\n')
velmenu_inp.close()

os.system('velmenu < velmenu.inp > /dev/null')

os.remove('velest.cmn')
os.remove('ttt.ttt')
os.remove('selstat.lis')
os.remove('invers.out')
os.remove('fin_hyp.cnv')
os.remove('data.nor')
os.remove('data.cnv')
os.remove('velmenu.inp')

d  = Read_Nordic(nordic_inp)
ce = Nordic2CNV(d)

# get refrence station


PRD(inp_file=nordic_inp,sta_inp='station.sta', falg_plot=False)


for _ in ['station.sta','model.mod','data.cnv']:
    
    if os.path.exists(os.path.join('..','velinp',_)):
   
        os.remove(os.path.join('..','velinp',_))

if not os.path.exists(os.path.join('..','velinp')):

    os.mkdir(os.path.join('..','velinp'))

move('station_ref.sta',os.path.join('..','velinp','station.sta'))
move('input.mod',os.path.join('..','velinp','initial.mod'))
move('data.cnv',os.path.join('..','velinp'))
os.remove('station.sta')


for i in ['gmap.cur.kml','hypmag.out','hypsum.out','print.out','hyp.out','tmp_sta','sta_cor.out']:
    
    if os.path.exists(i):
    
        os.remove(i)

print "\n+++ %d Events have been converted to CNV.\n"%(ce)
