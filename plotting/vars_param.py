# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
from  collections import defaultdict
import plot_settings as ps


"""
plot variables parameters
"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2017, UCAR/NOAA"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"



defs =  defaultdict(dict)

defs['elev']['var']   = 'zeta'
defs['elev']['vmin']  = -1.0
defs['elev']['vmax']  = 1.0
defs['elev']['label'] = 'Elev. [m]'
defs['elev']['format']= '%3.1g'
defs['elev']['cmap']  = plt.cm.jet

defs['elevdif']['var']   = 'zeta'
defs['elevdif']['vmin']  = -0.25
defs['elevdif']['vmax']  =  0.25
defs['elevdif']['label'] = 'Elev. Diff. [m]'
defs['elevdif']['format']= '%3.1g'
defs['elevdif']['cmap']  = ps.jetWoGn()

defs['dep']['var']   = 'h'
defs['dep']['vmin']  = 0.0
defs['dep']['vmax']  = 10 #6.0
defs['dep']['label'] = 'Depth [m]'
defs['dep']['format']= '%3.1g'
defs['dep']['cmap']  = plt.cm.jet_r

defs['temp']['var']   = 'temp'
defs['temp']['vmin']  = 22.0
defs['temp']['vmax']  = 25.0
defs['temp']['label'] = 'Temp. [.C]'
defs['temp']['format']= '%.1g'
defs['temp']['cmap']  = ps.jetWoGn()

defs['tempdif']['var']   = 'tempdif'
defs['tempdif']['vmin']  = -1.0
defs['tempdif']['vmax']  =  1.0
defs['tempdif']['label'] = 'Temp. diff. [.C]'
defs['tempdif']['format']= '%.1g'
defs['tempdif']['cmap']  = ps.jetWoGn()

defs['salt']['var']   = 'salt'
defs['salt']['vmin']  = 34.0
defs['salt']['vmax']  = 35.0
defs['salt']['label'] = 'Salt. [psu]'
defs['salt']['format']= '%.2g'
defs['salt']['cmap']  = ps.jetWoGn()

defs['saltdif']['var']   = 'saltdif'
defs['saltdif']['vmin']  = -1.0
defs['saltdif']['vmax']  =  1.0
defs['saltdif']['label'] = 'Salt diff. [psu]'
defs['saltdif']['format']= '%.2g'
defs['saltdif']['cmap']  = ps.jetWoGn()

defs['dens']['var']   = 'dens'
defs['dens']['vmin']  = 1021
defs['dens']['vmax']  = 1025
defs['dens']['label'] = 'density [kg/m \mathrm$^3$]'
defs['dens']['format']= '%.1g'
defs['dens']['cmap']  = ps.jetWoGn()

defs['vort']={}
defs['vort']['var']   = 'vort'
defs['vort']['vmin']  = -0.002
defs['vort']['vmax']  =  0.002
defs['vort']['label'] = 'Vorticity [s$^\mathrm{-1}$]'
defs['vort']['format']= '%.3g'
defs['vort']['cmap']  = ps.jetWoGn()

defs['div']['var']   = 'div'
defs['div']['vmin']  = -0.0015
defs['div']['vmax']  =  0.0015
defs['div']['label'] = 'Divergence [s$^\mathrm{-1}$]'
defs['div']['format']= '%.3g'
defs['div']['cmap']  = ps.jetWoGn()

defs['u']['var']   = 'u'
defs['u']['vmin']  = -0.5
defs['u']['vmax']  =  0.5
defs['u']['label'] = 'U [ms$^\mathrm{-1}$]'
defs['u']['format']= '%.2g'
defs['u']['cmap']  = ps.jetWoGn()

defs['v']['var']   = 'v'
defs['v']['vmin']  = -0.5
defs['v']['vmax']  =  0.5
defs['v']['label'] = 'V [ms$^\mathrm{-1}$]'
defs['v']['format']= '%.2g'
defs['v']['cmap']  = ps.jetWoGn()

defs['w']={}
defs['w']['var']   = 'w'
defs['w']['vmin']  = -0.5e-3
defs['w']['vmax']  =  0.5e-3
defs['w']['label'] = 'W [ms$^\mathrm{-1}$]'
defs['w']['format']= '%.4g'
defs['w']['cmap']  = ps.jetWoGn()

defs['uv']={}
defs['uv']['var']   = 'uv'
defs['uv']['vmin']  =  0.0
defs['uv']['vmax']  =  1.0
defs['uv']['label'] = 'Tot. Vel. [ms$^\mathrm{-1}$]'
defs['uv']['format']= '%.2g'
defs['uv']['cmap']  = ps.my_cmap

defs['discharge']={}
defs['discharge']['var']   = 'discharge'
defs['discharge']['vmin']  =  0.0
defs['discharge']['vmax']  =  30.0
defs['discharge']['label'] = 'Q [m$^\mathrm{3}$s$^\mathrm{-1}$]'
defs['discharge']['format']= '%.2g'
defs['discharge']['cmap']  = ps.my_cmap


defs['hs']={}
defs['hs']['var']   = 'Hwave'
defs['hs']['vmin']  =  0.
defs['hs']['vmax']  =  1.5 
defs['hs']['label'] = 'Hsig [m] '    
defs['hs']['format']= '%.1g'
defs['hs']['cmap']  = ps.my_cmap2

defs['wdir']={}
defs['wdir']['var']   = 'Dissip_wcap'
defs['wdir']['vmin']  =  0.
defs['wdir']['vmax']  =  360.0 
defs['wdir']['label'] = 'Wave mean Dir. [degN]'    
defs['wdir']['format']= '%.2g'
defs['wdir']['cmap']  = ps.my_cmap

defs['wcap']={}
defs['wcap']['var']   = 'Dwave'
defs['wcap']['vmin']  =  0.
defs['wcap']['vmax']  =  2.0 
defs['wcap']['label'] = 'Dissp. Wcap [Watm$^\mathrm{-2}$]'    
defs['wcap']['format']= '%.1g'
defs['wcap']['cmap']  = ps.my_cmap

defs['sin_wind']={}
defs['sin_wind']['var']   = 'win'
defs['sin_wind']['tvar']   = 'time'
defs['sin_wind']['vmin']  =  0.
defs['sin_wind']['vmax']  =  3e-5
defs['sin_wind']['label'] = 'Wind input [m$^ \mathrm{2}$s$^ \mathrm{-1}$] '
defs['sin_wind']['format']= '%.1g'
defs['sin_wind']['cmap']  = ps.my_cmap

defs['sdis_wave']={}
defs['sdis_wave']['var']   = 'dis'
defs['sdis_wave']['tvar']   = 'time'
defs['sdis_wave']['vmin']  =  0.
defs['sdis_wave']['vmax']  =  3e-5
defs['sdis_wave']['label'] = 'Surf. Dissip. [m$^ \mathrm{2}$s$^ \mathrm{-1}$] '
defs['sdis_wave']['format']= '%.1g'
defs['sdis_wave']['cmap']  = ps.my_cmap

defs['tke']={}
defs['tke']['var']   = 'tke'
defs['tke']['vmin']  = 0e-4
defs['tke']['vmax']  = 6e-4
defs['tke']['label'] = 'TKE [m$^ \mathrm{2}s$^ \mathrm{-2}]'
defs['tke']['format']= '%.1g'
defs['tke']['cmap']  = ps.jetWoGn()

defs['wlen']={}
defs['wlen']['var']   = 'Lwave'
defs['wlen']['vmin']  =  10.
defs['wlen']['vmax']  =  50. 
defs['wlen']['label'] = 'Wave length [m] '   
defs['wlen']['format']= '%.1g'
defs['wlen']['cmap']  = ps.my_cmap

defs['spm']={}
defs['spm']['var']   = 'SPM'
defs['spm']['vmin']  =  0.
defs['spm']['vmax']  =  2e-2
defs['spm']['label'] = 'SPM [kgm$^ \mathrm{-3}] '   
defs['spm']['format']= '%.1g'
defs['spm']['cmap']  = ps.my_cmap

defs['spmdif']={}
defs['spmdif']['var']   = 'SPM_diff'
defs['spmdif']['vmin']  =  -0.01
defs['spmdif']['vmax']  =   0.01
defs['spmdif']['label'] = 'SPM diff. [kgm$^ \mathrm{-3}] '   
defs['spmdif']['format']= '%.1g'
defs['spmdif']['cmap']  = ps.jetWoGn()


### Basic Defs
defs['wind']={}
defs['wind']['var']   = 'wind'
defs['wind']['tvar']   = 'time'
defs['wind']['vmin']  =  0.
defs['wind']['vmax']  =  10.
defs['wind']['label'] = 'Wind speed [ms$^ \mathrm{-1}$] '
defs['wind']['format']= '%.1g'
defs['wind']['cmap']  = ps.my_cmap


defs['windd']={}
defs['windd']['var']   = 'wdir.'
defs['windd']['tvar']   = 'time'
defs['windd']['vmin']  =  0.
defs['windd']['vmax']  =  360.
defs['windd']['label'] = 'Wind direction [deg] '
defs['windd']['format']= '%.1g'
defs['windd']['cmap']  = ps.my_cmap


#plot extent
defs['lim']['xmin']  = -100.0 
defs['lim']['xmax']  =  400.0
defs['lim']['ymin']  =  1000.0
defs['lim']['ymax']  =  2000.0
defs['lim']['zmin']  = -6.0
defs['lim']['zmax']  =  1.0


