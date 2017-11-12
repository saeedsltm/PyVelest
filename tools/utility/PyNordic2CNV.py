from PyNordicRW import Read_Nordic
from datetime import datetime as dt
from datetime import timedelta as td
from numpy import arange

"""
ChangeLog:

2017-Jun > Prevent repeating phases [The next same station and Phase ID will be ignored.]
2017-May > Now it's working for 19th and 20th century simultanously. 

"""

def Nordic2CNV(inp, out='data.cnv'):

    """
    inp: Nordic dictionary
    out: phase data in CNV format
    """

    fo = open(out,'w')
    ce = 0
    d  = inp 

    for evt in sorted(d.keys()):

        tmp_ph   = []
        a        = d[evt]
        ot_year  = a['HEADER']['L1']['Year']
        ot_month = a['HEADER']['L1']['Month']
        ot_day   = a['HEADER']['L1']['Day']
        ot_hour  = a['HEADER']['L1']['Hour']
        ot_min   = a['HEADER']['L1']['Min']
        ot_sec   = a['HEADER']['L1']['Sec']
        ot_msec  = int((ot_sec - int(ot_sec))*1e6)

        if ot_sec >= 60:

            ot_sec = ot_sec - 60
            ot_min+=1

        if ot_min >= 60:

            ot_min = ot_min - 60
            ot_hour+=1

        if ot_hour >= 24:

            ot_hour = ot_hour - 24
            ot_day+=1
                
        ot       = dt(ot_year,ot_month,ot_day,ot_hour,ot_min,int(ot_sec),ot_msec)
            
        lat      = a['HEADER']['L1']['Lat']
        lon      = a['HEADER']['L1']['Lon']
        dep      = a['HEADER']['L1']['Dep']
        mag      = a['HEADER']['L1']['Mag1']
        rms      = a['HEADER']['L1']['RMS']

        # check for gap
        
        if a['HEADER']['LE']:

            gap = a['HEADER']['LE']['GAP']

        else:

            gap = 999

        # check for mag
        
        if not str(mag).strip():

            mag = 0.0

        hdr = '%02d%02d%02d %02d%02d %5.2f %7.4fN  %7.4fE  %5.2f   %4.2f    %3s      %4.2f     '\
              %(int(str(ot_year)[-2:]), ot_month, ot_day, ot_hour, ot_min, ot_sec,lat, lon, dep,mag, gap, rms) 

        fo.write(hdr+'\n')
        
        body = []
        
        for i in a['PHASE']:

            tmp_flag = False
            otday    = ot_day

            #_P-PHASE

            for Pphas in a['PHASE'][i]['P']:

                if i+Pphas[0].upper() not in tmp_ph:

                    tmp_ph.append(i+Pphas[0].upper())

                else:

                    continue

                w = str(a['PHASE'][i]['P'][Pphas]['WI'])
                
                if not w.strip():

                    w = 0

                else:

                    w = int(w)
                
                try:

                    ar_hour = a['PHASE'][i]['P'][Pphas]['Hour']
                    ar_min  = a['PHASE'][i]['P'][Pphas]['Min']
                    ar_sec  = a['PHASE'][i]['P'][Pphas]['Sec']
                    ar_msec = int((ar_sec - int(ar_sec))*1e6)

                except (TypeError,ValueError):
                   
                    break

                if ar_sec >= 60:

                    ar_sec = ar_sec - 60
                    ar_min+=1

                if ar_min >= 60:

                    ar_min = ar_min - 60
                    ar_hour+=1

                if ar_hour >= 24:

                    ar_hour = ar_hour - 24
                    otday+=1
                        
                ar    = dt(ot_year,ot_month,otday,ar_hour,ar_min,int(ar_sec),ar_msec)  
                delta = ar - ot
                delta = delta.total_seconds()

                if delta > 200.0:

                    print '+++ Illogical phase time! (TT > 200sec), check > ',evt

                    break

                if delta > 0.0 and w < 4:
                    
                    tmp     = '%-4sP%1d%6.2f'%(i,w,delta)
                    body.append(tmp)

            otday    = ot_day
            
            #_S-PHASE

            for Sphas in a['PHASE'][i]['S']:

                if i+Sphas[0].upper() not in tmp_ph:

                    tmp_ph.append(i+Sphas[0].upper())

                else:

                    continue
                
                w = str(a['PHASE'][i]['S'][Sphas]['WI'])
                
                if not w.strip():

                    w = 0

                else:

                    w = int(w)
                    
                ar_hour = a['PHASE'][i]['S'][Sphas]['Hour']
                ar_min  = a['PHASE'][i]['S'][Sphas]['Min']
                ar_sec  = a['PHASE'][i]['S'][Sphas]['Sec']
                ar_msec = int((ar_sec - int(ar_sec))*1e6)
              
                if ar_sec >= 60:

                    ar_sec = ar_sec - 60
                    ar_min+=1

                if ar_min >= 60:

                    ar_min = ar_min - 60
                    ar_hour+=1

                if ar_hour >= 24:

                    ar_hour = ar_hour - 24
                    otday+=1

                ar    = dt(ot_year,ot_month,otday,ar_hour,ar_min,int(ar_sec),ar_msec)
                delta = ar - ot
                delta = delta.total_seconds()

                if delta > 200.0:

                    print '+++ Illogical phase time! (TT > 200sec), check > ',evt

                    break
                
                if delta > 0.0 and w < 4:
                    
                    tmp     = '%-4sS%1d%6.2f'%(i,w,delta)
                    body.append(tmp)

        complete_line   = len(body)/6
        incomplete_line = len(body)%6

        c = 0

        for i in range(complete_line):

            fo.write(''.join(body[0+c:6+c])+'\n')
            c+=6

            if i==complete_line-1:

                ce+=1

        if complete_line == 0 and incomplete_line:

            ce+=1

        if incomplete_line:

            fo.write(''.join(body[-incomplete_line:])+'\n')
            

        fo.write('\n')
    
    fo.close()

    print ''
    
    return ce


def ReadCNV(inpcnv='data.cnv'):

    cnv_dic = {}

    with open(inpcnv) as f:

        for l in f:

            if l.strip():

                chk_hdr = [h.isdigit() or h==' ' for h in l[0:5]]

                if all(chk_hdr):

                    yer  = int('%004d'%float(l[0:2]))
                    if 80 <= yer <= 99: yer+=1900
                    else: yer+=2000
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
    
def CNV2Nordic(inp='data.cnv',out='cnv.out'):

    """
    inp: phase data in CNV format.
    out: phase data in NORDIC format.
    """

    input_cnv  = ReadCNV(inp)   
    nordic_out = open(out,'w')

    for evt in input_cnv:

        l1  = ' %4d %02d%02d %02d%02d %4.1f L %7.3f %7.3f%5.1f                                    1\n'%(evt.year,evt.month,evt.day,
                                                                                                      evt.hour,evt.minute,evt.second+evt.microsecond*1e-6,
                                                                                                      input_cnv[evt]['lat'],input_cnv[evt]['lon'],input_cnv[evt]['dep'])

        nordic_out.write(l1)
        nordic_out.write(' STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU VELO AIN AR TRES W  DIS CAZ7\n')

        for i in range(len(input_cnv[evt]['arr']['sta'])):

            sta = input_cnv[evt]['arr']['sta'][i]
            pha = input_cnv[evt]['arr']['pha'][i]
            wte = input_cnv[evt]['arr']['wt'][i]
            art = evt + td(seconds=input_cnv[evt]['arr']['tim'][i])
            l4  = ' %-4s     %1s   %1s   %02d%02d %5.2f                                                   4\n'%(sta,pha,wte,
                                                                                                                art.hour,art.minute,art.second+art.microsecond*1e-6)
            nordic_out.write(l4)

        nordic_out.write('\n')

    nordic_out.close()

#___________________START
    
#d = Read_Nordic('select.out')
#Nordic2CNV(d)
#CNV2Nordic()
