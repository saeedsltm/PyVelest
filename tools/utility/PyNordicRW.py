from collections import OrderedDict

"""
A simple python script for reading SEISAN nordic
files.

Written By, S.SoltaniMoghadam__IIEES,Teh,2016

Notes:

1- It just read line types 1,4,E,F and the other lines
   will be ignored.
   
2- All variable names are consistance with SEISAN
nordic format. See SEISAN V10.0 manuall, appendix A
for more information.

Logs:

2017-Jun: Fix bugs for bad format in line-E, Line-4.
2017-Mar: Correct for station code when they are numbers.
2017-Mar: Program is renamed from "nordiceader" to "PyNordicRW".
2017-Mar: More than one phase [P/S] is allowded now for each station.
2017-Mar: Add a function to write NORDIC.
2017-Mar: Correct LINE dictionaries.

"""
#___________________ Define Nordic Type Lines as dictionaries:


line_1 = OrderedDict([
    ('F1',   {'col':(0,1),    'fmt':'%1s'    ,'v':' '  }), ('Year', {'col':(1,5),    'fmt':'%00004d','v':None }),
    ('F2',   {'col':(5,6),    'fmt':'%1s'    ,'v':' '  }), ('Month',{'col':(6,8),    'fmt':'%002d'  ,'v':None }),
    ('Day',  {'col':(8,10),   'fmt':'%002d'  ,'v':None }), ('FO',   {'col':(10,11),  'fmt':'%1s'    ,'v':' '  }),
    ('Hour', {'col':(11,13),  'fmt':'%002d'  ,'v':0    }), ('Min',  {'col':(13,15),  'fmt':'%002d'  ,'v':0    }),
    ('F3',   {'col':(15,16),  'fmt':'%1s'    ,'v':' '  }), ('Sec',  {'col':(16,20),  'fmt':'%4.1f'  ,'v':0    }),
    ('LMI',  {'col':(20,21),  'fmt':'%1s'    ,'v':' '  }), ('DI',   {'col':(21,22),  'fmt':'%1s'    ,'v':' '  }),
    ('EID',  {'col':(22,23),  'fmt':'%1s'    ,'v':' '  }), ('Lat',  {'col':(23,30),  'fmt':'%7.3f'  ,'v':None }),
    ('Lon',  {'col':(30,38),  'fmt':'%8.3f'  ,'v':None }), ('Dep',  {'col':(38,43),  'fmt':'%5.1f'  ,'v':None }),
    ('DepI', {'col':(43,44),  'fmt':'%1s'    ,'v':' '  }), ('LocI', {'col':(44,45),  'fmt':'%1s'    ,'v':' '  }),
    ('AGC',  {'col':(45,48),  'fmt':'%3s'    ,'v':'   '}), ('NumS', {'col':(48,51),  'fmt':'%0003d' ,'v':None }),
    ('RMS',  {'col':(51,55),  'fmt':'%4.1f'  ,'v':None }), ('Mag1', {'col':(55,59),  'fmt':'%4.1f'  ,'v':None }),
    ('TMag1',{'col':(59,60),  'fmt':'%1s'    ,'v':' '  }), ('MAGC1',{'col':(61,63),  'fmt':'%3s'    ,'v':'   '}),
    ('Mag2', {'col':(63,67),  'fmt':'%4.1f'  ,'v':None }), ('TMag2',{'col':(67,68),  'fmt':'%1s'    ,'v':' '  }),
    ('MAGC2',{'col':(68,71),  'fmt':'%3s'    ,'v':'   '}), ('MAG3', {'col':(71,75),  'fmt':'%4.1f'  ,'v':None }),
    ('TMag3',{'col':(75,76),  'fmt':'%1s'    ,'v':' '  }), ('MAGC3',{'col':(76,79),  'fmt':'%3s'    ,'v':'   '}),
    ('TL'  , {'col':(79,80),  'fmt':'%1s'    ,'v':' '  })])

line_E = OrderedDict([
    ('GAP', {'col':(5,8)  , 'fmt':'%3d'    ,'v':None}), ('F1' , {'col':(8,14) , 'fmt':'%6s'    ,'v':' ' }),
    ('ORT', {'col':(14,20), 'fmt':'%6.2f'  ,'v':None}), ('F2' , {'col':(20,24), 'fmt':'%4s'    ,'v':' ' }),
    ('LAE', {'col':(24,30), 'fmt':'%6.1f'  ,'v':None}), ('F3' , {'col':(30,32), 'fmt':'%2s'    ,'v':' ' }),
    ('LOE', {'col':(32,38), 'fmt':'%6.1f'  ,'v':None}), ('DEE', {'col':(38,44), 'fmt':'%5.1f'  ,'v':None}),
    ('CXY', {'col':(44,55), 'fmt':'%12.4e' ,'v':None}), ('CXZ', {'col':(55,67), 'fmt':'%12.4e' ,'v':None}),
    ('CYZ', {'col':(67,79), 'fmt':'%12.4e' ,'v':None})])

          
line_F = OrderedDict([
    ('STR',  {'col':(1,10) , 'fmt':'%10.4f' ,'v':None}), ('DIP' , {'col':(10,20), 'fmt':'%10.4f' ,'v':None}),
    ('RAK',  {'col':(20,30), 'fmt':'%10.4f' ,'v':None}), ('EHST', {'col':(30,35), 'fmt':'%5.1f'  ,'v':None}),
    ('EHD',  {'col':(35,40), 'fmt':'%5.1f'  ,'v':None}), ('EHR' , {'col':(40,45), 'fmt':'%5.1f'  ,'v':None}),
    ('EF1',  {'col':(45,50), 'fmt':'%5.1f'  ,'v':None}), ('EF2' , {'col':(50,55), 'fmt':'%5.1f'  ,'v':None}),
    ('AMR',  {'col':(55,60), 'fmt':'%5.1f'  ,'v':None}), ('NBP' , {'col':(60,62), 'fmt':'%2d'    ,'v':None}),
    ('NBA' , {'col':(62,64), 'fmt':'%2d'    ,'v':None}), ('F1'  , {'col':(64,66), 'fmt':'%2s'    ,'v':' ' }),
    ('AGC',  {'col':(66,69), 'fmt':'%3s'    ,'v':' ' }), ('PRO' , {'col':(70,77),  'fmt':'%7s'   ,'v':' ' }),
    ('QUA',  {'col':(77,80), 'fmt':'%3s'    ,'v':' ' })])

line_6 = {'fn':   {'col':(1,79)    ,'fmt':'%79s'  ,'v':None}}

line_4 = OrderedDict([
    ('F1'  , {'col':(0,1)  , 'fmt':'%1s'  ,'v':' ' }), ('STA' , {'col':(1,6)  , 'fmt':'%-5s'  ,'v':' ' }),
    ('INST', {'col':(6,7)  , 'fmt':'%1s'  ,'v':' ' }), ('COMP', {'col':(7,8)  , 'fmt':'%1s'   ,'v':' ' }),
    ('W'   , {'col':(8,9)  , 'fmt':'%1s'  ,'v':' ' }), ('QLTI', {'col':(9,10) , 'fmt':'%1s'   ,'v':' ' }),
    ('PH'  , {'col':(10,14), 'fmt':'%-4s' ,'v':' ' }), ('WI'  , {'col':(14,15), 'fmt':'%1s'   ,'v':' ' }),
    ('F2'  , {'col':(15,16), 'fmt':'%1s'  ,'v':' ' }), ('POL' , {'col':(16,17), 'fmt':'%1s'   ,'v':' ' }),
    ('F3'  , {'col':(17,18), 'fmt':'%1s'  ,'v':' ' }), ('Hour', {'col':(18,20), 'fmt':'%002d' ,'v':None}),
    ('Min' , {'col':(20,22), 'fmt':'%002d','v':None}), ('Sec' , {'col':(22,28), 'fmt':'%6.2f' ,'v':None}),
    ('F4'  , {'col':(28,29), 'fmt':'%1s'  ,'v':' ' }), ('DTN' , {'col':(29,33), 'fmt':'%4d'   ,'v':None}),
    ('AMP' , {'col':(33,40), 'fmt':'%7.1f','v':None}), ('F5'  , {'col':(40,41), 'fmt':'%1s'   ,'v':' ' }),
    ('PER' , {'col':(41,45), 'fmt':'%4.0f','v':None}), ('F6'  , {'col':(45,46), 'fmt':'%1s'   ,'v':' ' }),
    ('DOA' , {'col':(46,51), 'fmt':'%5.0f','v':None}), ('F7'  , {'col':(51,52), 'fmt':'%1s'   ,'v':' ' }),
    ('PVEL', {'col':(52,56), 'fmt':'%4.0f','v':None}), ('ANGI', {'col':(56,60), 'fmt':'%4.0f' ,'v':None}),
    ('ARES', {'col':(60,63), 'fmt':'%3d'  ,'v':None}), ('TRES', {'col':(63,68), 'fmt':'%5.1f' ,'v':None}),
    ('WI2' , {'col':(68,70), 'fmt':'%2d'  ,'v':None}), ('DIS' , {'col':(70,75), 'fmt':'%5.0f' ,'v':None}),
    ('F8'  , {'col':(75,76), 'fmt':'%1s'  ,'v':' ' }), ('AZ'  , {'col':(76,79), 'fmt':'%3d'   ,'v':None}),
    ('TL'  , {'col':(79,80), 'fmt':'%1s'  ,'v':' ' })])


Flag_1 = True
Flag_4 = True
Flag_6 = False
Flag_E = True
Flag_F = False

def Read_Nordic(inp='select.out'):

    """
    Input: 'select.out' or any nordic file.
    Output: a dictionary containing the followings keys/values:

    out_dic:
           |
           |
            Event_ID_dic:[ex: 20130105123412.7]
                        |
                        |
                         Station_dic: [ex: AAA]
                                    |
                                    |
                                     Phase_ID: [ex: P]
                                             |
                                             |
                                              key:val [ex: 'Lat':54.43, ...]
    """

    nor_dic = OrderedDict()
    
    with open(inp) as f:

        flag_l4 = False
        flag_l1 = True

        #__________ Sec 01, Read Line Type 1:

        for l in f:

            if l[79:80] == '1' and flag_l1:

                flag_l4 = False
                flag_l1 = False

                EID = '%4d%002d%002d%002d%002d%04.1f'%(float(l[1:5]), float(l[6:8]), float(l[8:10]),
                                                  float(l[11:13]), float(l[13:15]), float(l[16:20]))
                nor_dic[EID] = OrderedDict()
                nor_dic[EID]['HEADER'] = OrderedDict()
                nor_dic[EID]['HEADER']['L1'] = OrderedDict()
                nor_dic[EID]['HEADER']['LE'] = OrderedDict()
                nor_dic[EID]['HEADER']['LF'] = OrderedDict()
                nor_dic[EID]['PHASE'] = OrderedDict()
                
                for i in line_1:

                    a,b = line_1[i]['col']

                    if l[a:b].strip():

                        u = l[a:b]

                        if '.' in u:

                            u = float(u.strip())

                        elif u.strip().isdigit():

                            u = int(u.strip())

                        else:

                            u = u.strip()
                    else:

                        fmt = line_1[i]['fmt']
                        fmt = fmt.split('%')[1].split(fmt[-1])[0].split('.')[0]
                        u   = '%'+fmt+'s'
                        u   = u%(' ')

                    nor_dic[EID]['HEADER']['L1'][i] = u

        #__________ Sec 02, Read Line Type E:

            if l[79:80] == 'E':

                for k in line_E:

                    a,b = line_E[k]['col']

                    if l[a:b].strip():

                        z = l[a:b]

                        if '.' in z:

                            try:
                                z = float(z.strip())

                            except ValueError:

                                print '+++ Skipping Header Value! Check Line-E at event => %s'%(EID)

                                pass

                        elif z.strip().isdigit():

                            z = int(z.strip())

                        else:

                            z = z.strip()
                    else:

                        fmt = line_E[k]['fmt']
                        fmt = fmt.split('%')[1].split(fmt[-1])[0].split('.')[0]
                        z   = '%'+fmt+'s'
                        z   = z%(' ')

                    nor_dic[EID]['HEADER']['LE'][k] = z
                            
        #__________ Sec 03, Read Line Type F:

            if Flag_F and l[79:80] == 'F':

                for p in ['FOCMEC','FPFIT','PINV']:

                    if p in l[70:77]:

                        nor_dic[EID]['HEADER']['LF'][p] = OrderedDict()

                        for k in line_F:

                            a,b = line_F[k]['col']

                            if l[a:b].strip():

                                z = l[a:b]

                                if '.' in z:

                                    z = float(z.strip())

                                elif z.strip().isdigit():

                                    z = int(z.strip())

                                else:

                                    z = z.strip()
                            else:

                                fmt = line_F[k]['fmt']
                                fmt = fmt.split('%')[1].split(fmt[-1])[0].split('.')[0]
                                z   = '%'+fmt+'s'
                                z   = z%(' ')
                                
                            nor_dic[EID]['HEADER']['LF'][p][k] = z

                        nor_dic[EID]['HEADER']['LF'][p]['PRO'] = p                   

        #__________ Sec 04, Read Line Type 4:

            if l[1:5] == 'STAT':

                flag_l4 = True

            if not l[1:5].strip():

                flag_l4 == False

            if flag_l4 and l[1:5] != 'STAT':

                flag_l1 = True

                s = l[1:5].strip()

                if s and s not in nor_dic[EID]['PHASE'].keys():

                    nor_dic[EID]['PHASE'][s] = {'P':OrderedDict(),'S':OrderedDict()}

                for j in line_4:

                    a,b = line_4[j]['col']

                    if l[a:b].strip():

                        v = l[a:b]

                        if '.' in v:

                            try:

                                v = float(v.strip())
                                
                            except ValueError:

                                print '+++ Skipping Phase! Check Line-4 at event => %s'%(EID)

                                pass
                    
                        elif (a,b) != (1,6) and v.strip().isdigit():

                            v = int(v.strip())

                        else:

                            v = v.strip()

                    else:

                        fmt = line_4[j]['fmt']
                        fmt = fmt.split('%')[1].split(fmt[-1])[0].split('.')[0]
                        v   = '%'+fmt+'s'
                        v   = v%(' ')

                    if 'P' == l[10:11].upper():

                        if l[10:14].strip() not in nor_dic[EID]['PHASE'][s]['P']:

                            nor_dic[EID]['PHASE'][s]['P'][l[10:14].strip()] = OrderedDict()

                        nor_dic[EID]['PHASE'][s]['P'][l[10:14].strip()][j] = v

                    if 'S' == l[10:11].upper():

                        if l[10:14].strip() not in nor_dic[EID]['PHASE'][s]['S']:

                            nor_dic[EID]['PHASE'][s]['S'][l[10:14].strip()] = OrderedDict()

                        nor_dic[EID]['PHASE'][s]['S'][l[10:14].strip()][j] = v
                        
##    print '+++ Read %d events successfully.'%(len(nor_dic))

    return nor_dic


def Write_Nordic(inp=None, output='pyselect.out'):

    output = open(output,'w')

    for evt in inp:

        header = inp[evt]['HEADER']
        phase  = inp[evt]['PHASE']

        #_____LINE 1

        for item in header['L1']:

            # FILL HEADER (LINE_1) WITH WITHE SPACE IF NO VALUE IS REPORTED.

            try:

                output.write(line_1[item]['fmt']%(header['L1'][item]))

            except TypeError:

                output.write(header['L1'][item])

        output.write('\n')

        #_____LINE E

        if header['LE']: output.write(' GAP=')

        for item in header['LE']:

            # FILL HEADER (LINE_E) WITH WITHE SPACE IF NO VALUE IS REPORTED.

            try:

                output.write(line_E[item]['fmt']%(header['LE'][item]))

            except TypeError:

                output.write(header['LE'][item])

        if header['LE']: output.write('E\n')

        #_____LINE F

        for fm in header['LF']:

            # FILL HEADER (LINE_F) WITH WITHE SPACE IF NO VALUE IS REPORTED.

            for item in header['LF'][fm]:

                try:

                    output.write(line_F[item]['fmt']%(header['LF'][fm][item]))

                except TypeError:

                    output.write(header['LF'][fm][item])

        if header['LF']: output.write('F\n')

        #_____LINE 7

        output.write(' STAT SP IPHASW D HRMN SECON CODA AMPLIT PERI AZIMU VELO AIN AR TRES W  DIS CAZ7\n')

        #_____LINE 4

        for sta in phase:

            for pha in phase[sta]:

                for ps in phase[sta][pha]:

                    for item in phase[sta][pha][ps]:

                        try:

                            output.write(line_4[item]['fmt']%(phase[sta][pha][ps][item]))

                        except TypeError:

                            fmt = line_4[item]['fmt']
                            fmt = '%'+str(int(float(fmt.split('%')[-1][:-1])))+'s'
                            
                            output.write(fmt%(' '))

                    output.write('\n')

        output.write('\n')
        
    output.close()

# EXAMPLE

#d = Read_Nordic('../initial.out')            
#Write_Nordic(inp=d, output='select.out')


