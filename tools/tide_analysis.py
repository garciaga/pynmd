__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Portland State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"

import  netCDF4
import  netcdftime
from    datetime import datetime,timedelta
from    dateutil.parser import *
import  scipy.io as sio
from    collections  import defaultdict
import  numpy as np
import  os,sys
import  pandas as pd
import  cPickle as pickle
import  string
import  matplotlib.pyplot as plt
import  pynmd.plotting.plot_settings as ps
import  glob

def projection_nri ():
    from mpl_toolkits.basemap import Basemap
    proj    = Basemap(resolution=None,projection='npstere',lon_0=-45,boundinglat=70)
    return proj
#
def lonlat2xy(proj,lon,lat,inverse=False,lon_orig= -77.3382,lat_orig=34.5279 ):
    xu_orig,yu_orig=proj(lon_orig, lat_orig)
    if not inverse:
        xu, yu = proj(lon, lat)
        x = xu - xu_orig
        y = yu - yu_orig
        return x,y
    else:
        x = lon + xu_orig
        y = lat + yu_orig        
        lonu, latu = proj(x,y,inverse=True)
        return lonu, latu       
#
def get_depth(x,y,all=False):
    bat_in   = 'data/nri_regional_grd_roms_v2.nc'
    tnc    = netCDF4.Dataset(bat_in)
    tvars  = tnc.variables
    xb     = tnc.variables['x_rho'][:]
    yb     = tnc.variables['y_rho'][:]
    bat    = tnc.variables['h'][:]
    bat    = np.ma.masked_where(bat==-10.0,bat)
    if False:
        dep    = calc.griddata(xb[~bat.mask], yb[~bat.mask],bat[~bat.mask] , x, y,extrap=True)
    else:
        import octant.csa as csa
        csa_interp = csa.CSA(xb[~bat.mask], yb[~bat.mask],bat[~bat.mask])
        dep = csa_interp(x,y)
    
    if all:
        return dep,xb,yb,bat
    else:
        return dep

def read_Elgar_pressure_data(flow_data):
    print '   > Elgar pressure data;'
    dir1 = data_dir + '/elgar/*_filter_abs.nc'
    flist    = glob.glob(dir1)
    flist.sort()

    for filename in flist[:]:
        ncf   = netCDF4.Dataset(filename,'r')
        ncvar = ncf.variables
        time       = ncvar['time']
        utime      = netcdftime.utime(time.units)
        dates      = utime.num2date(time[:])
        sta_name   = filename.split('/')[-1][:3]
        if 'p' in sta_name or 'q' in sta_name:
            print '    > read  > Station name: ', sta_name
    
            water_depth = ncvar['water_depth'][:]
    
            data  = pd.DataFrame(data = water_depth, columns = ['water_depth'], index = dates)    
            data  = data.dropna()
            data  = data.resample('H')  #hourly mean
            
            date_tmp =  ps.datetime64todatetime(data.index.values)
            datenum  = [ps.datetime2datenum(datei) for datei in  date_tmp]

            
            flow_data[sta_name]['name_long'] = 'Bottom pressure '+sta_name
            flow_data[sta_name]['lat']       = ncvar['lat'][:]
            flow_data[sta_name]['lon']       = ncvar['lon'][:]
            flow_data[sta_name]['elev']      = (data['water_depth'] - data['water_depth'].mean()).values
            flow_data[sta_name]['date']      = date_tmp
            flow_data[sta_name]['datenum_mat']   = datenum

#
def read_wl_noaa_station(flow_data,inp_dir,sta_name):
    print  inp_dir
    dir1  = data_dir + inp_dir + '/'
    finfo = open(dir1 + 'info')
    for line in finfo:
        words = string.split(line)
        if 'lat ' in line: lat = float(words[1])
        if 'lon ' in line: lon = float(words[1])
    finfo.close()
    
    filename2 = dir1 + 'data.csv'
    print 'Read observation at: ', filename2
    
    fp   = open(filename2, "r")
    line = ''
    line = fp.readline()
    
    obs = []
    obs_dates = []
    for line in fp:
        words = string.split(line,',')
        #Date Time, Water Level, Sigma, I, L 
        try:
           obs_dates.append(parse(words[0]))
        except:
           break
        wl = float(words[1])            # m
        obs.append(wl)
    
    fp.close()
    
    obs       = np.array(obs)
    obs_dates = np.array(obs_dates)
    
    print '    > read  > Station name: ', sta_name
    
    water_depth = obs
    data  = pd.DataFrame(data = water_depth, columns = ['water_depth'], index = obs_dates)    
    data  = data.dropna()
    data  = data.resample('H')  #hourly mean
    
    date_tmp =  ps.datetime64todatetime(data.index.values)
    datenum  = [ps.datetime2datenum(datei) for datei in  date_tmp]
    
    flow_data[sta_name]['name_long'] = inp_dir
    flow_data[sta_name]['lat']       = lat
    flow_data[sta_name]['lon']       = lon
    flow_data[sta_name]['elev']      = (data['water_depth'] - data['water_depth'].mean()).values
    flow_data[sta_name]['date']      = date_tmp
    flow_data[sta_name]['datenum_mat']   = datenum

#
def read_ncom(flow_data):
    print '  > NCOM;'
    filename  = data_dir + '/ncom/ncom_glb_regp01_2012-3-4-5.nc'

    ncf   = netCDF4.Dataset(filename,'r')
    ncvar = ncf.variables
    time       = ncvar['time']
    utime      = netcdftime.utime(time.units)
    dates      = utime.num2date(time[:])
    sta_name   = filename.split('/')[-1][:3]
    
    lona,lata = np.meshgrid(ncvar['lon'][:],ncvar['lat'][:])
    dist = ((lona - lon_orig)**2 + (lata - lat_orig)**2)
    [i,j] = np.where(dist==dist.min())
    
    elev  = ncvar['surf_el'][:,i,j].flatten()
    data  = pd.DataFrame(data = elev, columns = ['elev'], index = dates)    
    data  = data.dropna()
    data  = data.resample('H')  #hourly mean
    
    date_tmp =  ps.datetime64todatetime(data.index.values)
    datenum  = [ps.datetime2datenum(datei) for datei in  date_tmp]
    
    flow_data[sta_name]['name_long'] = 'Elevation from NCOM model '
    flow_data[sta_name]['lat']       = lata[i,j].item()
    flow_data[sta_name]['lon']       = lona[i,j].item()
    flow_data[sta_name]['elev']      = data['elev'].values
    flow_data[sta_name]['date']      = date_tmp
    flow_data[sta_name]['datenum_mat']   = datenum
#
def read_roms_zeta(flow_data,name,roms_hisfile):
    print '  >  ROMS ZETA  > in   read_roms_zeta();'
    ncvar   = netCDF4.Dataset(roms_hisfile,'r').variables
    dates   = netCDF4.num2date(times=ncvar['ocean_time'][:],units=ncvar['ocean_time'].units,calendar='standard')
    x       = ncvar['x_rho'][:]
    y       = ncvar['y_rho'][:]
    elev    = ncvar['zeta'][:]
    nj,ni   = elev[0,:].shape
    mask    = elev[0,:].mask
    
    for i in range(ni):
        for j in range(nj):
            if ~mask[j,i]:
                sta_name   = name+'__nj_0'+str(j)+'__ni_0'+str(i)+'_'
                #print sta_name
                data  = pd.DataFrame(data = elev[:,j,i], columns = ['elev'], index = dates)    
                data  = data.dropna()
                data  = data.resample('H')  #hourly mean
                
                date_tmp =   ps.datetime64todatetime(data.index.values)
                datenum  =  [ps.datetime2datenum(datei) for datei in  date_tmp]
                
                flow_data[sta_name]['case']      = name
                flow_data[sta_name]['name_long'] = 'Elevation from ROMS model '+sta_name.replace('_',' ')
                flow_data[sta_name]['i']         = i
                flow_data[sta_name]['j']         = j
                flow_data[sta_name]['y']         = y[j,i].item()
                flow_data[sta_name]['x']         = x[j,i].item()                
                flow_data[sta_name]['elev']      = data['elev'].values
                flow_data[sta_name]['date']      = date_tmp
                flow_data[sta_name]['datenum_mat']   = datenum
            #
    proj = projection_nri ()
    #add xy
    for sta_name in flow_data.keys():
        flow_data[sta_name]['lon'],flow_data[sta_name]['lat'] = lonlat2xy(proj,flow_data[sta_name]['x'],flow_data[sta_name]['y'],inverse=True)
    

def mat2py_datenum (mtime):
    """
    A serial date number represents a calendar date as the number of
    days that has passed since a fixed base date. In MATLAB, serial
    date number 1 is January 1, 0000.
    """
    date=[]
    for it in range(len(mtime)):
        date.append(datetime.fromordinal(int(mtime[it])) + timedelta(days=mtime[it]%1) - timedelta(days = 366))
    return pylab.array(date)


def datetime64todatetime(dt):
    tmp=[]
    for it in range(len(dt)):
        tmp.append(pd.Timestamp(dt[it]).to_pydatetime())
    return pylab.array(tmp)


def datetime2datenum(dt):
    """
    PY-datetime 2 Matlab datenum (by saeed)
    No idea why adding 367 need to figure out! maybe not right!!!    :-(
    """
    delta = dt - datetime (1,1,1)
    return delta.days + 367.0 + delta.seconds / 86400.0


def datetime2matlabdn(dt):
   """
   from:
   http://stackoverflow.com/questions/8776414/python-datetime-to-matlab-datenum
   """
   ord = dt.toordinal()
   mdn = dt + timedelta(days = 366)
   frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   return mdn.toordinal() + frac


def read_noaa_station(fname):
    print 'Read observation for: ', fname
    #
    fp   = open(fname, "r")
    line = ''
    line = fp.readline()
    #
    obs = []
    obs_dates = []
    for line in fp:
        words = string.split(line)
        #Date Time, Water Level, Sigma, I, L 
        try:
           
           tstr =  words[1]+'-'+words[2]+'-'+words[3]+' '+words[4]+':'+words[5]
           obs_dates.append( parse(tstr))
        except:
           break
        wl = float(words[-1])            # m
        obs.append(wl)
    #
    fp.close()
    #
    obs       = np.array(obs)
    obs_dates = np.array(obs_dates)
    #
    data  = pd.DataFrame(data = obs, columns = ['wl'], index = obs_dates)    
    data  = data.dropna()
    #data  = data.resample('H')  #hourly mean

    date_tmp =  datetime64todatetime(data.index.values)
    return dict(dates=date_tmp, wl = data.wl.values)



#
def do_r_t_tide_analysis(flow_data,constits,out_dir):
    """
    only works for hourly data
    use resample function in dataframe to resample it to 'H'
    
    TODO: remove hourly constraint
    """
    r_t_tide_out_dir = out_dir + '/r_t_tide_out_txt_files/' 
    os.system('mkdir -p '+ r_t_tide_out_dir)
     
    #prepare pickle name !
    pick_dir = r_t_tide_out_dir + '/r_t_tide_pickle/' 
    os.system('mkdir -p '+ pick_dir)
    pick_name = pick_dir + 'rt_tide_data.pickle'   
   
    if (not os.path.exists(pick_name)):
        import matlab.engine
        eng = matlab.engine.start_matlab()
        eng.addpath(r'/data01/01-projects/01-passaic/03-modeling/01-delft3d/02-pass-ideal/pycodes/post/06-r_t_tide/mcodes/r_t_tide',nargout=0);
        #
        for sta_name in flow_data.keys():
            print sta_name

            t_sta  = np.array(flow_data[sta_name]['datenum_mat'])
            h_sta  = np.array(flow_data[sta_name]['elev'])
            [rk]   = np.where( ~np.isnan(h_sta) )
            t_tmp  = t_sta [rk].tolist()
            h_tmp  = h_sta [rk].tolist()
            out_file = r_t_tide_out_dir +flow_data[sta_name]['name_long'].replace(' ','_')+'.rtt'
            
            #run matlab r_t_tide
            
            try:
                lat =  flow_data[sta_name]['lat'].item()
            except:
                lat = flow_data[sta_name]['lat']

            s = eng.r_t_tide_saeed    (matlab.double(t_tmp),matlab.double(h_tmp),'latitude',
                lat ,'nodalcorrflag','true','greenwichcorrflag','true','method','cauchy',
                'output',out_file,'rayleigh',0.9 );
            
            #vectors to pick the requested constituents
            sta_cnam   = []
            sta_amp    = []
            sta_amp_ci = []
            sta_pha    = []
            sta_pha_ci = []
            
            #read t_tide_text file
            f_rtt = open(out_file)
            nam    = []
            amp    = []
            amp_ci = []
            pha    = []
            pha_ci = []
            start_read = False
            for line in f_rtt.readlines():
                if 'pha_err' in line: 
                    start_read = True
                    continue
                if start_read:
                    tmp = line.split()
                    nam.append(line[1:5].rstrip())
                    amp.append(   float(tmp[2]))
                    amp_ci.append(float(tmp[3]))
                    pha.append(   float(tmp[4]))
                    pha_ci.append(float(tmp[5]))
            f_rtt.close()
            
            nam    = np.array(nam)
            amp    = np.array(amp)
            amp_ci = np.array(amp_ci)
            pha    = np.array(pha)
            pha_ci = np.array(pha_ci)
            
            ##    
            for const in  constits:
                sta_cnam.append   (const)
                if const in nam:
                    [ic] = np.where(const == nam)
                    sta_amp.append   (amp   [ic])
                    sta_amp_ci.append(amp_ci[ic])
                    sta_pha.append   (pha   [ic])
                    sta_pha_ci.append(pha_ci[ic])
                else:
                    cmask = np.array([9.999e-12])
                    sta_amp.append   (cmask)
                    sta_amp_ci.append(cmask)
                    sta_pha.append   (cmask)
                    sta_pha_ci.append(cmask)                    #if len(const) > 2:
            
            sta_cnam    = np.array(sta_cnam).squeeze()  
            sta_amp     = np.array(sta_amp).squeeze() 
            sta_amp_ci  = np.array(sta_amp_ci).squeeze() 
            sta_pha     = np.array(sta_pha).squeeze() 
            sta_pha_ci  = np.array(sta_pha_ci).squeeze() 
            sta_all = np.array(zip(sta_cnam,sta_amp,sta_amp_ci,sta_pha,sta_pha_ci))                    
            sta_df  = pd.DataFrame(data=sta_all, columns=['Const','amp','amp_err','pha','pha_err'])
            
            sta_df['amp'        ] = sta_df['amp'    ].convert_objects(convert_numeric=True)
            sta_df['amp_err'    ] = sta_df['amp_err'].convert_objects(convert_numeric=True)
            sta_df['pha'        ] = sta_df['pha'    ].convert_objects(convert_numeric=True)
            sta_df['pha_err'    ] = sta_df['pha_err'].convert_objects(convert_numeric=True)
            
            flow_data[sta_name]['r_t_tide'] = sta_df
            
            #delete elev data
            del flow_data[sta_name]['elev']
        
        print ' > Write pickle > ', pick_name
        pickle.dump( flow_data, open(pick_name , "wb" ) )
    else:
        print 'Read pickle > ', pick_name
        flow_data = pickle.load( open( pick_name , "r" ) )        
    
    return r_t_tide_out_dir        
    
def do_tappy_tide_analysis(flow_data,constits,out_dir):
    from tappy_local import tappy   
    
    tappy_tide_out_dir = out_dir + '/tappy_tide_out_txt_files/' 
    os.system('mkdir -p '+ tappy_tide_out_dir)    

    #prepare pickle name !
    pick_dir = tappy_tide_out_dir + '/tappy_tide_pickle/' 
    os.system('mkdir -p '+ pick_dir)
    pick_name = pick_dir + 'tappy_tide_data.pickle'
    
    if (not os.path.exists(pick_name)):
        for sta_name in flow_data.keys():
            dates      = np.array(flow_data[sta_name]['date'])
            elev       = np.array(flow_data[sta_name]['elev'])
            lon        = np.array(flow_data[sta_name]['lon'])
            lat        = np.array(flow_data[sta_name]['lat'])
            
            out_file = tappy_tide_out_dir + flow_data[sta_name]['name_long'].replace(' ','_')+'.tappy'
            
            
            ### Saeed tries to understand!  from here
            data_filename    ='test'
            def_filename     = None
            config           = None
            quiet            = False
            debug            = False
            outputts         = False
            outputxml        = ''
            ephemeris        = False
            rayleigh         = 0.9
            print_vau_table  = False
            missing_data     = 'ignore'
            #missing_data    = 'fill'
            linear_trend     = False
            remove_extreme   = False
            zero_ts          = None
            filter           = None
            pad_filters      = None
            include_inferred = True
            xmlname          = flow_data[sta_name]['name_long']
            xmlcountry       = 'US'
            xmllatitude      = lat
            xmllongitude     = lon
            xmltimezone      = '0000'
            xmlcomments      = 'No comment'
            xmlunits         = 'm or ms-1'
            xmldecimalplaces = None
            
            ############## model
            x = tappy(
                outputts  = outputts,
                outputxml = 'model.xml',
                quiet     = quiet,
                debug     = debug,
                ephemeris = ephemeris,
                rayleigh  = rayleigh,
                print_vau_table = print_vau_table,
                missing_data = missing_data,
                linear_trend = linear_trend,
                remove_extreme = remove_extreme,
                zero_ts = zero_ts,
                filter  = filter,
                pad_filters = pad_filters,
                include_inferred = include_inferred,
                )
            
            x.dates      = dates
            x.elevation  = elev
            package      = x.astronomic(x.dates)
            (x.zeta, x.nu, x.nup, x.nupp, x.kap_p, x.ii, x.R, x.Q, x.T, x.jd, x.s, x.h, x.N, x.p, x.p1) = package
            ray = 1.0
            (x.speed_dict, x.key_list) = x.which_constituents(len(x.dates),package,rayleigh_comp = ray)
            
            x.constituents()
            x.print_con()
            x.print_con_file(filedat = out_file, lon = lon, lat = lat)
            only_const = False
            
            if not only_const:
                if x.missing_data == 'fill':
                        x.dates_filled, x.elevation_filled = x.missing(x.missing_data, x.dates, x.elevation)
                        x.write_file( x.dates_filled,
                                        x.elevation_filled,
                                        fname='outts_filled.dat')
                
                x.filter = 'usgs'
                #x.filter='doodson'
                
                if x.filter:
                        for item in x.filter.split(','):
                            if item in ['mstha', 'wavelet', 'cd', 'boxcar', 'usgs', 'doodson', 'lecolazet1', 'kalman', 'transform']:# 'lecolazet', 'godin', 'sfa']:
                                filtered_dates, result = x.filters(item, x.dates, x.elevation)
                                x.write_file(filtered_dates, result, fname='outts_filtered_%s.dat' % (item,))
                                x_dates_filter= filtered_dates
                                x_eleva_filter= result              
                        (x.speed_dict, x.key_list) = x.which_constituents(len(x.dates),package,rayleigh_comp = ray)
                
            #vectors to pick the requested constituents
            sta_cnam   = []
            sta_amp    = []
            sta_pha    = []
            
            #read t_tide_text file
            f_rtt = open(out_file)
            nam    = []
            amp    = []
            pha    = []
            start_read = False
            for line in f_rtt.readlines():
                if '=        =====' in line: 
                    start_read = True
                    continue
                
                if len(line) < 3:
                    break
            
                if start_read:
                    tmp = line.split()
                    nam.append(line[9:14].strip())
                    amp.append(   float(tmp[2]))
                    pha.append(   float(tmp[3]))
            
            f_rtt.close()
            
            nam    = np.array(nam)
            amp    = np.array(amp)
            pha    = np.array(pha)
            
            ##    
            for const in  constits:
                sta_cnam.append   (const)
                if const in nam:
                    [ic] = np.where(const == nam)
                    sta_amp.append   (amp   [ic].item())
                    sta_pha.append   (pha   [ic].item())
                else:
                    cmask = np.array([9.999e-12])
                    sta_amp.append   (cmask)
                    sta_pha.append   (cmask)
            
            sta_cnam = np.array(sta_cnam).flatten()
            sta_amp  = np.array(sta_amp ).flatten()
            sta_pha  = np.array(sta_pha ).flatten()  
            
            sta_all  = np.array(zip(sta_cnam,sta_amp,sta_pha))                    
            sta_df   = pd.DataFrame(data = sta_all, columns=['Const','amp','pha'])
            
            sta_df['amp'] = sta_df['amp'].convert_objects(convert_numeric=True)
            sta_df['pha'] = sta_df['pha'].convert_objects(convert_numeric=True)
            
            flow_data[sta_name]['tappy_tide'] = sta_df
            print ' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '
        

        print ' > Write pickle > ', pick_name
        pickle.dump( flow_data, open(pick_name , "wb" ) )
        #
    else:
        print 'Read pickle > ', pick_name
        flow_data = pickle.load( open( pick_name , "r" ) )        
    
    return tappy_tide_out_dir        
    
    
 
