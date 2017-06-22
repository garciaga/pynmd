__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Portland State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"

import numpy as np
import os,sys

def write_init(fname,data):
    """
    Write Delft3D init txt file
    - data    : a dict contan required data (elev,u,v,temp,salt)
        u   (nz,ny,nx)
        v   (nz,ny,nx)
        temp(nz,ny,nx)
        salt(nz,ny,nx)

    - filename  
    
    
    #We need all this fields in one text file
    #1) water level
    #2) u-velocity (for each (k)-layer)
    #3) v-velocity (for each (k)-layer)
    #4) salinity (for each (k)-layer, if selected)
    #5) temperature (for each (k)-layer)
    #6) constituent 1(for each (k)-layer)
    #7) constituent 2 (for each (k)-layer)
    #8) constituent n (for each (k)-layer)
    """
    
    os.system('rm '+ fname)
    # read main inputs
    elev = data['elev'][:]
    uu   = data['uu']  [:]
    vv   = data['vv']  [:]
    salt = data['salt'][:]
    temp = data['temp'][:]
    
    if 'sed' in data.keys():
       sed  = data['sed'][:]
    
    #
    nlayer =   temp.shape[0]       
    ## Write
    fmt = '%.4g'
    f   = open(fname,'a')
    
    #1)water level
    np.savetxt(f,elev.T , fmt=fmt)
    
    #2) u-velocity (for each (k)-layer)
    for il in range (nlayer):
        np.savetxt(f,uu[il].T, fmt=fmt)
    
    #3) v-velocity (for each (k)-layer)
    for il in range (nlayer):
        np.savetxt(f,vv[il].T, fmt=fmt)
    
    #4) salinity (for each (k)-layer, if selected)
    for il in range (nlayer):
        np.savetxt(f,salt[il].T, fmt=fmt)
    
    #5) Temp (for each (k)-layer, if selected)
    for il in range (nlayer):
        np.savetxt(f,temp[il].T, fmt=fmt)
   
    if 'sed' in data.keys():
    #6) one sed component  (for each (k)-layer, if selected)
        if len(sed.shape) > 3 : 
            if len(sed.shape) > 3 : 
                for ised in range (len(sed)): #sys.exit('ERR : only one sed comp is implemented')
                    for il in range (nlayer):
                        np.savetxt(f,sed[ised,il].T, fmt=fmt)
        else:        
            for il in range (nlayer):
                np.savetxt(f,sed[il].T, fmt=fmt)
    
    f.close()
    
    

def write_rgh(fname,data):
    """
    Write Delft3D roughness txt file
    - data    : a dict contan required data (elev,u,v,temp,salt)
        rfgx (ny,nx)
        rfgy

    - filename  
    
    """
    os.system('rm '+ fname)
    # read main inputs
    rghx = data['rghx'][:]
    rghy = data['rghy'][:]
    #
    ## Write
    fmt = '%.4g'
    f   = open(fname,'a')
    
    #1)roughness for ux and uy
    np.savetxt(f,rghx.T , fmt=fmt)
    np.savetxt(f,rghy.T , fmt=fmt)
   
   
    f.close()
    
def write_dep(fname,data):
    """
    Write Delft3D dep txt file
    - data    : a dict contan  dep (ny,nx)

    - filename  
    
    """
    os.system('rm '+ fname)
    # read main inputs
    dep = data['dep'][:]
    #
    ## Write
    fmt = '%.4g'
    f   = open(fname,'a')
    
    #1)roughness for ux and uy
    np.savetxt(f,dep.T , fmt=fmt)
  
    f.close()        

def write_table(fname,data):
    """
    To write one table block of text boundary into .bcc or .bct file:
    
    In some cases you need to read in the file into GUI and save it again
    to be able to read it by Delft3D properly.
    
    Linear : time interpolation  
    Linear : Depth interpolation (e.g salinity)
    """
    fmt = "%.7e"
    f   = open(fname,"a")
    
    params_header = []
    params_header.append( "table-name          '"  +  data["table_name"]    +"'\n") 
    params_header.append( "contents            '"  +  data["contents"]      +"'\n") 
    params_header.append( "location            '"  +  data["location"]      +"'\n") 
    params_header.append( "time-function       '"  +  data["time-function"] +"'\n") 
    params_header.append( "reference-time       "  +  data["time-ref"]      +" \n") 
    params_header.append( "time-unit           '"  +  data["time-unit"]     +"'\n") 
    params_header.append( "interpolation       '"  +  data["interpolation"] +"'\n") 
    
    main_param      = data["main_param"]
    main_param_unit = data["main_param_unit"]
    params_header.append(     "parameter           'time                '  unit '[min]' \n")
    if "linear" in data["contents"].lower():
        params_header.append( "parameter           '" +  main_param +  "   end A surface'       unit '" + main_param_unit+"' \n")
        params_header.append( "parameter           '" +  main_param +  "   end A bed    '       unit '" + main_param_unit+"' \n")
        params_header.append( "parameter           '" +  main_param +  "   end B surface'       unit '" + main_param_unit+"' \n")
        params_header.append( "parameter           '" +  main_param +  "   end B bed    '       unit '" + main_param_unit+"' \n")
    elif "uniform" in data["contents"].lower():     
        params_header.append( "parameter           '" +  main_param +  "  end A uniform'        unit '" + main_param_unit+"' \n")
        params_header.append( "parameter           '" +  main_param +  "  end B uniform'        unit '" + main_param_unit+"' \n")
    else:
        sys.exit(data["contents"]+' > contents < method not yet implemented !!')
    
    data = data["data"]
    nl,_ = data.shape
    params_header.append( "records-in-table "  +  str(int(nl)) +"\n") 
    
    for il in range(len(params_header)):
        f.write(params_header[il])
    
    np.savetxt(f,data , fmt=fmt)
    f.close()
    #####


def discharge_salt_sed():
    ref_time  = datetime.datetime(2010,01,01)
    ref_units = "minutes since "+ref_time.isoformat()
    ut       = netcdftime.utime(ref_units.replace("T"," "))    
    
    mins  = np.arange(0,3.6820800e+006 * 2, 3.6820800e+006)
    dates = ut.num2date(mins)
    #RIVER    
    #write discharge file
    data = {}
    fname                   = "ExpPass1.bct"    
    os.system('rm -rf '+fname)

    data["table_name"]      = "Boundary Section : 2"
    data["contents"]        = "Uniform             "
    data["location"]        = "River               "
    data["time-ref"]        = ref_time.isoformat()[:10].replace("-","")
    data["time-function"]   = "non-equidistant"
    data["time-unit"]       = ut.units
    data["interpolation"]   = "linear"
    data["main_param"]      = "total discharge (t)"
    data["main_param_unit"] = "[m3/s]"
    discharge = np.ones((len(mins))) * -17.0
    timevec   = ut.date2num(dates)
    data["data"]            = np.array(zip(timevec,discharge,np.ones_like(discharge)*999.999))
    write_table(fname=fname,data=data)

    #write BCC (e.g. salinity) file
    #SeaSide
    #salinity
    data = {}
    fname                   = "ExpPass1.bcc"    
    os.system('rm -rf '+fname)
    data["table_name"]      = "Boundary Section : 1"
    data["contents"]        = "Linear             "
    data["location"]        = "SeaSide            "
    data["time-function"]   = "non-equidistant"
    data["time-ref"]        = ref_time.isoformat()[:10].replace("-","")
    data["time-unit"]       = ut.units
    data["interpolation"]   = "linear"
    #
    data["main_param"]      = "Salinity"
    data["main_param_unit"] = "[ppt]"
    #end A
    end_a_srf = np.ones((len(mins))) * 17.0
    end_a_bot = np.ones((len(mins))) * 25.0
    #end B
    end_b_srf = np.ones((len(mins))) * 17.0
    end_b_bot = np.ones((len(mins))) * 25.0
    #
    timevec   = ut.date2num(dates)
    data["data"] = np.array(zip(timevec,end_a_srf,end_a_bot,end_b_srf,end_b_bot))
    write_table(fname=fname,data=data)
    #
    #Sediment
    data = {}
    fname                   = "ExpPass1.bcc"    
    data["table_name"]      = "Boundary Section : 1"
    data["contents"]        = "Uniform            "
    data["location"]        = "SeaSide            "
    data["time-function"]   = "non-equidistant"
    data["time-ref"]        = ref_time.isoformat()[:10].replace("-","")
    data["time-unit"]       = ut.units
    data["interpolation"]   = "linear"
    #
    data["main_param"]      = "Sediment_maximum"
    data["main_param_unit"] = "[kg/m3]"
    end_a = np.ones((len(mins))) * 0.0  #end A
    end_b = np.ones((len(mins))) * 0.0  #end B 
    timevec   = ut.date2num(dates)
    data["data"]            = np.array(zip(timevec,end_a,end_b))
    write_table(fname=fname,data=data)
    #river side
    #SeaSide
    #salinity
    data = {}
    data = {}
    fname                   = "ExpPass1.bcc"    
    data["table_name"]      = "Boundary Section : 2"
    data["contents"]        = "Uniform            "
    data["location"]        = "River              "
    data["time-function"]   = "non-equidistant"
    data["time-ref"]        = ref_time.isoformat()[:10].replace("-","")
    data["time-unit"]       = ut.units
    data["interpolation"]   = "linear"
    #
    data["main_param"]      = "Salinity"
    data["main_param_unit"] = "[ppt]"
    end_a = np.ones((len(mins))) * 0.0  #end A
    end_b = np.ones((len(mins))) * 0.0  #end B 
    timevec   = ut.date2num(dates)
    data["data"]            = np.array(zip(timevec,end_a,end_b))
    write_table(fname=fname,data=data)
    #
    #Sediment
    data = {}
    fname                   = "ExpPass1.bcc"    
    data["table_name"]      = "Boundary Section : 2"
    data["contents"]        = "Uniform            "
    data["location"]        = "River              "
    data["time-function"]   = "non-equidistant"
    data["time-ref"]        = ref_time.isoformat()[:10].replace("-","")
    data["time-unit"]       = ut.units
    data["interpolation"]   = "linear"
    #
    data["main_param"]      = "Sediment_maximum"
    data["main_param_unit"] = "[kg/m3]"
    end_a = np.ones((len(mins))) * 0.02  #end A
    end_b = np.ones((len(mins))) * 0.02  #end B 
    timevec   = ut.date2num(dates)
    data["data"]            = np.array(zip(timevec,end_a,end_b))
    write_table(fname=fname,data=data)


def discharge_salt():
    ref_time  = datetime.datetime(2010,01,01)
    ref_units = "minutes since "+ref_time.isoformat()
    ut       = netcdftime.utime(ref_units.replace("T"," "))    
    
    mins  = np.arange(0,3.6820800e+006, 3.6820800e+006 // 60)
    dates = ut.num2date(mins)
    #RIVER    
    #write discharge file
    data = {}
    fname                   = "ExpPass1.bct"    
    os.system('rm -rf '+fname)

    data["table_name"]      = "Boundary Section : 2"
    data["contents"]        = "Uniform             "
    data["location"]        = "River               "
    data["time-ref"]        = ref_time.isoformat()[:10].replace("-","")
    data["time-function"]   = "non-equidistant"
    data["time-unit"]       = ut.units
    data["interpolation"]   = "linear"
    data["main_param"]      = "total discharge (t)"
    data["main_param_unit"] = "[m3/s]"
    discharge = np.ones((len(mins))) * -17.0
    timevec   = ut.date2num(dates)
    data["data"]            = np.array(zip(timevec,discharge,np.ones_like(discharge)*999.999))
    write_table(fname=fname,data=data)

    #write BCC (e.g. salinity) file
    #SeaSide
    #salinity
    data = {}
    fname                   = "ExpPass1.bcc"    
    os.system('rm -rf '+fname)
    data["table_name"]      = "Boundary Section : 1"
    data["contents"]        = "Linear             "
    data["location"]        = "SeaSide            "
    data["time-function"]   = "non-equidistant"
    data["time-ref"]        = ref_time.isoformat()[:10].replace("-","")
    data["time-unit"]       = ut.units
    data["interpolation"]   = "linear"
    #
    data["main_param"]      = "Salinity"
    data["main_param_unit"] = "[ppt]"
    #end A
    end_a_srf = np.ones((len(mins))) * 17.0
    end_a_bot = np.ones((len(mins))) * 25.0
    #end B
    end_b_srf = np.ones((len(mins))) * 17.0
    end_b_bot = np.ones((len(mins))) * 25.0
    #
    timevec   = ut.date2num(dates)
    data["data"] = np.array(zip(timevec,end_a_srf,end_a_bot,end_b_srf,end_b_bot))
    write_table(fname=fname,data=data)
    #
    #river side
    #salinity
    data = {}
    data = {}
    fname                   = "ExpPass1.bcc"    
    data["table_name"]      = "Boundary Section : 2"
    data["contents"]        = "Uniform            "
    data["location"]        = "River              "
    data["time-function"]   = "non-equidistant"
    data["time-ref"]        = ref_time.isoformat()[:10].replace("-","")
    data["time-unit"]       = ut.units
    data["interpolation"]   = "linear"
    #
    data["main_param"]      = "Salinity"
    data["main_param_unit"] = "[ppt]"
    end_a = np.ones((len(mins))) * 0.0  #end A
    end_b = np.ones((len(mins))) * 0.0  #end B 
    timevec   = ut.date2num(dates)
    data["data"]            = np.array(zip(timevec,end_a,end_b))
    write_table(fname=fname,data=data)
    #


if __name__ == "__main__":
    # execute only if run as a script
    discharge_salt()

# # 
# #     BCT    
# table-name           'Boundary Section : 2'
# contents             'Uniform             '
# location             'River               '
# time-function        'non-equidistant'
# reference-time       20100401
# time-unit            'minutes'
# interpolation        'linear'
# parameter            'time                '                     unit '[min]'
# parameter            'total discharge (t)  end A'               unit '[m3/s]'
# parameter            'total discharge (t)  end B'               unit '[m3/s]'
# records-in-table     2
#  0.0000000e+000 -1.7000000e+001  9.9999900e+002
#  4.1760000e+004 -1.7000000e+001  9.9999900e+002
# 
# 
# 
# 
# 
# 
