from LatLon import lat_lon as ll

#___________________READ NORDIC STATION0.HYP FILE INTO A DICTIONARY.

def Read_Nordic_Sta(inp='STATION0.HYP'):

    """

    Read NORDIC STATION0.HYP and parse it
    into a dictionary.

    """

    sta_dic = {'RESET':[],
               'STA':{},
               'VELO':{'Vp':[],'Dep':[]},
               'CNTL':{'ST_D':None, 'X0':None,'X1':None, 'VpVs':None,'NumSD':None,'STR':None,'INC':None}
               }

    with open(inp) as f:

        flag = 0

        for l in f:

            if not l.strip():

                flag+=1

            if l.startswith('RESET'):

                sta_dic['RESET'].append(l.strip())

            elif len(l)>22 and l[13] in ['N','S'] and l[22] in ['E','W']:

                sta_nam = l[:6].strip()
                sta_lat = ll.Latitude(degree=float(l[6:8]),minute=float(l[8:13])).decimal_degree
                sta_lon = ll.Longitude(degree=float(l[14:17]),minute=float(l[17:22])).decimal_degree
                sta_elv = float(l[23:27])
                pcr     = None
                scr     = None

                if len(l) >= 33:

                    pcr = float(l[27:33])

                if len(l) >= 40:

                    scr = float(l[34:40])
                    
                if sta_nam not in sta_dic['STA']:

                    sta_dic['STA'][sta_nam] = {'LAT':sta_lat,
                                               'LON':sta_lon,
                                               'ELV':sta_elv,
                                               'PCR':pcr,
                                               'SCR':scr
                                               }
            if flag == 2 and l.split():

                vd = l.split()

                sta_dic['VELO']['Vp'].append(float(vd[0]))
                sta_dic['VELO']['Dep'].append(float(vd[1]))

            if flag == 3 and l.strip():

                sta_dic['CNTL']['ST_D']=float(l[:5])
                sta_dic['CNTL']['X0']=float(l[5:10])
                sta_dic['CNTL']['X1']=float(l[10:15])
                sta_dic['CNTL']['VpVs']=float(l[15:20])

                try:

                    sta_dic['CNTL']['NumSD']=float(l[20:25])
                    sta_dic['CNTL']['STR']=float(l[25:30])
                    sta_dic['CNTL']['INC']=float(l[30:35])

                except (IndexError,ValueError):

                    pass

                break

    return sta_dic

#___________________WRITE NORDIC STATION0.HYP DICTIONARY INTO A NORIC FORMAT.

def Write_Nordic_Sta(inp_dic, output='STATION0_NEW.HYP'):

    output = open(output,'w')

    for reset in inp_dic['RESET']:

        output.write(reset+'\n')

    output.write('\n')
    
    for sta in sorted(inp_dic['STA']):

        lon = ll.Longitude(inp_dic['STA'][sta]['LON'])
        lat = ll.Latitude(inp_dic['STA'][sta]['LAT'])
        elv = inp_dic['STA'][sta]['ELV']
        pcr = inp_dic['STA'][sta]['PCR']
        scr = inp_dic['STA'][sta]['SCR']

        if pcr != None and scr != None:

            output.write('  %-4s%02d%05.2f%1s%3d%05.2f%1s%00004d %+5.2f  %+5.2f\n'%(sta,
                                                                                    lat.degree,lat.decimal_minute,lat.get_hemisphere(),
                                                                                    lon.degree,lon.decimal_minute,lon.get_hemisphere(),
                                                                                    elv,pcr,scr))
        elif pcr != None and scr == None:

            output.write('  %-4s%02d%05.2f%1s%3d%05.2f%1s%00004d %+5.2f\n'%(sta,
                                                                            lat.degree,lat.decimal_minute,lat.get_hemisphere(),
                                                                            lon.degree,lon.decimal_minute,lon.get_hemisphere(),
                                                                            elv,pcr))            

        else:

            output.write('  %-4s%02d%05.2f%1s%3d%05.2f%1s%00004d\n'%(sta,
                                                                       lat.degree,lat.decimal_minute,lat.get_hemisphere(),
                                                                       lon.degree,lon.decimal_minute,lon.get_hemisphere(),
                                                                       elv))            
    output.write('\n')

    for vel,dep in zip(sorted(inp_dic['VELO']['Vp']),sorted(inp_dic['VELO']['Dep'])):

        output.write('%7.3f%7.3f\n'%(vel,dep))

    output.write('\n')

    st_d = inp_dic['CNTL']['ST_D']
    x_in = inp_dic['CNTL']['X0']
    x_fi = inp_dic['CNTL']['X1']
    vpvs = inp_dic['CNTL']['VpVs']
    nmsd = inp_dic['CNTL']['NumSD']
    strd = inp_dic['CNTL']['STR']
    inc  = inp_dic['CNTL']['INC']

    if nmsd != None and strd != None and inc != None:

        output.write('%4.0f.%4.0f.%4.0f.%5.2f%4.0f.%4.0f.%4.0f.\n'%(st_d,x_in,x_fi,vpvs,nmsd,strd,inc))

    else:

        output.write('%4.0f.%4.0f.%4.0f.%5.2f\n'%(st_d,x_in,x_fi,vpvs))

    output.write('NEW')
    output.close()

#___________________START
    
##inp = Read_Nordic_Sta()
##Write_Nordic_Sta(inp)        
