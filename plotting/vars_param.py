# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
from  collections import defaultdict
import colormaps as maps


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
defs['elev']['vmin']  = -1
defs['elev']['vmax']  =  1
defs['elev']['label'] = 'Elev. [m]'
defs['elev']['format']= '%3.1g'
defs['elev']['cmap']  = maps.jetWoGn()

defs['elevdif']['var']   = 'zeta'
defs['elevdif']['vmin']  = -0.25
defs['elevdif']['vmax']  =  0.25
defs['elevdif']['label'] = 'Elev. Diff. [m]'
defs['elevdif']['format']= '%3.1g'
defs['elevdif']['cmap']  = maps.jetWoGn()

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
defs['temp']['cmap']  = maps.jetWoGn()

defs['tempdif']['var']   = 'tempdif'
defs['tempdif']['vmin']  = -1.0
defs['tempdif']['vmax']  =  1.0
defs['tempdif']['label'] = 'Temp. diff. [.C]'
defs['tempdif']['format']= '%.1g'
defs['tempdif']['cmap']  = maps.jetWoGn()

defs['salt']['var']   = 'salt'
defs['salt']['vmin']  = 34.0
defs['salt']['vmax']  = 35.0
defs['salt']['label'] = 'Salt. [psu]'
defs['salt']['format']= '%.3g'
defs['salt']['cmap']  = maps.jetWoGn()

defs['saltdif']['var']   = 'saltdif'
defs['saltdif']['vmin']  = -1.0
defs['saltdif']['vmax']  =  1.0
defs['saltdif']['label'] = 'Salt diff. [psu]'
defs['saltdif']['format']= '%.2g'
defs['saltdif']['cmap']  = maps.jetWoGn()

defs['dens']['var']   = 'dens'
defs['dens']['vmin']  = 1021
defs['dens']['vmax']  = 1025
defs['dens']['label'] = 'density [kg/m \mathrm$^3$]'
defs['dens']['format']= '%.1g'
defs['dens']['cmap']  = maps.jetWoGn()

defs['vort']={}
defs['vort']['var']   = 'vort'
defs['vort']['vmin']  = -0.002
defs['vort']['vmax']  =  0.002
defs['vort']['label'] = 'Vorticity [s$^\mathrm{-1}$]'
defs['vort']['format']= '%.3g'
defs['vort']['cmap']  = maps.jetWoGn()

defs['div']['var']   = 'div'
defs['div']['vmin']  = -0.0015
defs['div']['vmax']  =  0.0015
defs['div']['label'] = 'Divergence [s$^\mathrm{-1}$]'
defs['div']['format']= '%.3g'
defs['div']['cmap']  = maps.jetWoGn()

defs['u']['var']   = 'u'
defs['u']['vmin']  = -0.5
defs['u']['vmax']  =  0.5
defs['u']['label'] = 'U [ms$^\mathrm{-1}$]'
defs['u']['format']= '%.2g'
defs['u']['cmap']  = maps.jetWoGn()

defs['v']['var']   = 'v'
defs['v']['vmin']  = -0.5
defs['v']['vmax']  =  0.5
defs['v']['label'] = 'V [ms$^\mathrm{-1}$]'
defs['v']['format']= '%.2g'
defs['v']['cmap']  = maps.jetWoGn()

defs['w']={}
defs['w']['var']   = 'w'
defs['w']['vmin']  = -0.5e-3
defs['w']['vmax']  =  0.5e-3
defs['w']['label'] = 'W [ms$^\mathrm{-1}$]'
defs['w']['format']= '%.4g'
defs['w']['cmap']  = maps.jetWoGn()

defs['x']={}
defs['x']['var']   = 'x'
defs['x']['vmin']  = -5e3
defs['x']['vmax']  =  5e3
defs['x']['label'] = 'X [m]'
defs['x']['format']= '%.4g'
defs['x']['cmap']  = maps.jetWoGn()

defs['y']={}
defs['y']['var']   = 'y'
defs['y']['vmin']  = -5e3
defs['y']['vmax']  =  5e3
defs['y']['label'] = 'Y [m]'
defs['y']['format']= '%.4g'
defs['y']['cmap']  = maps.jetWoGn()

defs['z']={}
defs['z']['var']   = 'z'
defs['z']['vmin']  = -5e3
defs['z']['vmax']  =  5e3
defs['z']['label'] = 'Z [m]'
defs['z']['format']= '%.4g'
defs['z']['cmap']  = maps.jetWoGn()

defs['surface']={}
defs['surface']['var']   = 'surface'
defs['surface']['vmin']  = -5e3
defs['surface']['vmax']  =  5e3
defs['surface']['label'] = 'Surf. [m$^\mathrm{2}$]'
defs['surface']['format']= '%.4g'
defs['surface']['cmap']  = maps.jetWoGn()

defs['volume']={}
defs['volume']['var']   = 'volume'
defs['volume']['vmin']  = -5e3
defs['volume']['vmax']  =  5e3
defs['volume']['label'] = 'Vol. [m$^\mathrm{3}$]'
defs['volume']['format']= '%.4g'
defs['volume']['cmap']  = maps.jetWoGn()

defs['uv']={}
defs['uv']['var']   = 'uv'
defs['uv']['vmin']  =  0.0
defs['uv']['vmax']  =  1.0
defs['uv']['label'] = 'Tot. Vel. [m s$^\mathrm{-1}$]'
defs['uv']['format']= '%.2g'
defs['uv']['cmap']  = maps.jetMinWi

defs['discharge']={}
defs['discharge']['var']   = 'discharge'
defs['discharge']['vmin']  =  0.0
defs['discharge']['vmax']  =  30.0
defs['discharge']['label'] = 'Q [m$^\mathrm{3}$ s$^\mathrm{-1}$]'
defs['discharge']['format']= '%.2g'
defs['discharge']['cmap']  = maps.jetMinWi


defs['hs']={}
defs['hs']['var']    = 'Hwave'
defs['hs']['vmin']   =  0.
defs['hs']['vmax']   =  1.5 
defs['hs']['label']  = 'Hsig [m] '    
defs['hs']['format'] = '%.1g'
defs['hs']['cmap']   = maps.jetMinWi2

defs['wdir']={}
defs['wdir']['var']    = 'Dwave'
defs['wdir']['vmin']   =  0.
defs['wdir']['vmax']   =  360.0 
defs['wdir']['label']  = 'Wave mean Dir. [degN]'    
defs['wdir']['format'] = '%.2g'
defs['wdir']['cmap']   = maps.jetMinWi2

defs['wcap']={}
defs['wcap']['var']    = 'Dissip_wcap'
defs['wcap']['vmin']   =  0.
defs['wcap']['vmax']   =  1.e-9 
defs['wcap']['label']  = 'Dissp. Wcap [Wat m$^\mathrm{-2}$]'    
defs['wcap']['format'] = '%.1g'
defs['wcap']['cmap']   = maps.jetMinWi2

defs['sin_wind']={}
defs['sin_wind']['var']    = 'win'
defs['sin_wind']['tvar']   = 'time'
defs['sin_wind']['vmin']   =  0.
defs['sin_wind']['vmax']   =  3e-5
defs['sin_wind']['label']  = 'Wind input [m$^ \mathrm{2}$ s$^ \mathrm{-1}$] '
defs['sin_wind']['format'] = '%.1g'
#defs['sin_wind']['cmap']   = maps.jetMinWi

defs['sdis_wave']={}
defs['sdis_wave']['var']   = 'dis'
defs['sdis_wave']['tvar']   = 'time'
defs['sdis_wave']['vmin']  =  0.
defs['sdis_wave']['vmax']  =  3e-5
defs['sdis_wave']['label'] = 'Surf. Dissip. [m$^ \mathrm{2}$ s$^ \mathrm{-1}$] '
defs['sdis_wave']['format']= '%.1g'
defs['sdis_wave']['cmap']  = maps.jetMinWi

defs['tke']={}
defs['tke']['var']   = 'tke'
defs['tke']['vmin']  = 0e-4
defs['tke']['vmax']  = 6e-4
defs['tke']['label'] = 'TKE [m$^ \mathrm{2}s$^ \mathrm{-2}]'
defs['tke']['format']= '%.1g'
defs['tke']['cmap']  = maps.jetWoGn()

defs['wlen']={}
defs['wlen']['var']   = 'Lwave'
defs['wlen']['vmin']  =  10.
defs['wlen']['vmax']  =  50. 
defs['wlen']['label'] = 'Wave length [m] '   
defs['wlen']['format']= '%.1g'
defs['wlen']['cmap']  = maps.jetMinWi

defs['spm']={}
defs['spm']['var']   = 'SPM'
defs['spm']['vmin']  =  0.
defs['spm']['vmax']  =  2e-2
defs['spm']['label'] = 'SPM [kgm$^ \mathrm{-3}] '   
defs['spm']['format']= '%.1g'
defs['spm']['cmap']  = maps.jetMinWi

defs['spmdif']={}
defs['spmdif']['var']   = 'SPM_diff'
defs['spmdif']['vmin']  =  -0.01
defs['spmdif']['vmax']  =   0.01
defs['spmdif']['label'] = 'SPM diff. [kg m$^ \mathrm{-3}] '   
defs['spmdif']['format']= '%.1g'
defs['spmdif']['cmap']  = maps.jetWoGn()


### Basic Defs
defs['wind']={}
defs['wind']['var']   = 'wind'
defs['wind']['tvar']   = 'time'
defs['wind']['vmin']  =  0.
defs['wind']['vmax']  =  10.
defs['wind']['label'] = 'Wind speed [m s$^ \mathrm{-1}$] '
defs['wind']['format']= '%.1g'
defs['wind']['cmap']  = maps.jetMinWi


defs['windd']={}
defs['windd']['var']   = 'wdir.'
defs['windd']['tvar']   = 'time'
defs['windd']['vmin']  =  0.
defs['windd']['vmax']  =  360.
defs['windd']['label'] = 'Wind direction [deg] '
defs['windd']['format']= '%.1g'
defs['windd']['cmap']  = maps.jetMinWi


### ROMS momentum diagnostics
defs['cor']={}
defs['cor']['var']   = 'cor'
defs['cor']['vmin']  =  -2e-5
defs['cor']['vmax']  =   2e-5
defs['cor']['label'] = 'Coriolis [m s$^ \mathrm{-2}$] '   
defs['cor']['format']= '%.2g'
defs['cor']['cmap']  = maps.jetWoGn()

defs['fsco']={}
defs['fsco']['var']   = 'fsco'
defs['fsco']['vmin']  = -1e-5
defs['fsco']['vmax']  =  1e-5
defs['fsco']['label'] ='horizontal Coriolis-stokes [m s$^ \mathrm{-2}$] '   
defs['fsco']['format']= '%.2g'
defs['fsco']['cmap']  = maps.jetWoGn()

defs['vadv']={}
defs['vadv']['var']   = 'vadv'
defs['vadv']['vmin']  =  -0.0002
defs['vadv']['vmax']  =   0.0002
defs['vadv']['label'] ='vertical advection [m s$^ \mathrm{-2}$] '   
defs['vadv']['format']= '%.2g'
defs['vadv']['cmap']  = maps.jetWoGn()

defs['hadv']={}
defs['hadv']['var']   = 'hadv'
defs['hadv']['vmin']  =  -0.0002
defs['hadv']['vmax']  =   0.0002
defs['hadv']['label'] ='horizontal advection [m s$^ \mathrm{-2}$] '   
defs['hadv']['format']= '%.2g'
defs['hadv']['cmap']  = maps.jetWoGn()

defs['xadv']={}
defs['xadv']['var']   = 'xadv'
defs['xadv']['vmin']  =  -0.01
defs['xadv']['vmax']  =   0.01
defs['xadv']['label'] ='horizontal XI-advection [m s$^ \mathrm{-2}$] '   
defs['xadv']['format']= '%.2g'
defs['xadv']['cmap']  = maps.jetWoGn()

defs['yadv']={}
defs['yadv']['var']   = 'yadv'
defs['yadv']['vmin']  =  -0.01
defs['yadv']['vmax']  =   0.01
defs['yadv']['label'] ='horizontal ETA-advection [m s$^ \mathrm{-2}$] '   
defs['yadv']['format']= '%.2g'
defs['yadv']['cmap']  = maps.jetWoGn()

defs['vjvf']={}
defs['vjvf']['var']   = 'vjvf'
defs['vjvf']['vmin']  =  -5e-6
defs['vjvf']['vmax']  =   5e-6
defs['vjvf']['label'] ='vertical J vortex force [m s$^ \mathrm{-2}$] '   
defs['vjvf']['format']= '%.2g'
defs['vjvf']['cmap']  = maps.jetWoGn()

defs['hjvf']={}
defs['hjvf']['var']   = 'hjvf'
defs['hjvf']['vmin']  =  -1e-5
defs['hjvf']['vmax']  =   1e-5
defs['hjvf']['label'] ='J vortex force [m s$^ \mathrm{-2}$] '   
defs['hjvf']['format']= '%.2g'
defs['hjvf']['cmap']  = maps.jetWoGn()

defs['kvrf']={}
defs['kvrf']['var']   = 'kvrf'
defs['kvrf']['vmin']  =  -1e-5
defs['kvrf']['vmax']  =   1e-5
defs['kvrf']['label'] ='K vortex force [m s$^ \mathrm{-2}$] '   
defs['kvrf']['format']= '%.2g'
defs['kvrf']['cmap']  = maps.jetWoGn()

defs['wbrk']={}
defs['wbrk']['var']   = 'wbrk'
defs['wbrk']['vmin']  =  -1e-8
defs['wbrk']['vmax']  =   1e-8
defs['wbrk']['label'] ='Wave breaking [m s$^ \mathrm{-2}$] '   
defs['wbrk']['format']= '%.2g'
defs['wbrk']['cmap']  = maps.jetWoGn()

defs['wrol']={}
defs['wrol']['var']   = 'wrol'
defs['wrol']['vmin']  =  -0.01
defs['wrol']['vmax']  =   0.01
defs['wrol']['label'] ='wave roller acceleration [m s$^ \mathrm{-2}$] '   
defs['wrol']['format']= '%.2g'
defs['wrol']['cmap']  = maps.jetWoGn()

defs['prsgrd']={}
defs['prsgrd']['var']   = 'prsgrd'
defs['prsgrd']['vmin']  =  -0.0002
defs['prsgrd']['vmax']  =   0.0002
defs['prsgrd']['label'] ='pressure gradient [m s$^ \mathrm{-2}$] '   
defs['prsgrd']['format']= '%.2g'
defs['prsgrd']['cmap']  = maps.jetWoGn()

defs['vvisc']={}
defs['vvisc']['var']   = 'vvisc'
defs['vvisc']['vmin']  =   0.0
defs['vvisc']['vmax']  =   1e-5
defs['vvisc']['label'] ='vertical viscosity [m s$^ \mathrm{-2}$] '   
defs['vvisc']['format']= '%.2g'
defs['vvisc']['cmap']  = maps.jetWoGn()

defs['accel']={}
defs['accel']['var']   = 'accel'
defs['accel']['vmin']  =  -0.00005
defs['accel']['vmax']  =   0.00005
defs['accel']['label'] ='acceleration [m s$^ \mathrm{-2}$] '   
defs['accel']['format']= '%.2g'
defs['accel']['cmap']  = maps.jetWoGn()

defs['salt_hadv']={}
defs['salt_hadv']['var']   = 'salt_hadv'
defs['salt_hadv']['vmin']  =  -0.05
defs['salt_hadv']['vmax']  =   0.05
defs['salt_hadv']['label'] ='Salt horizontal advection [s$^ \mathrm{-1}$] '   
defs['salt_hadv']['format']= '%.2g'
defs['salt_hadv']['cmap']  = maps.jetWoGn()

defs['salt_vadv']={}
defs['salt_vadv']['var']   = 'salt_vadv'
defs['salt_vadv']['vmin']  =  -0.05
defs['salt_vadv']['vmax']  =   0.05
defs['salt_vadv']['label'] ='Salt vertical advection [s$^ \mathrm{-1}$] '   
defs['salt_vadv']['format']= '%.2g'
defs['salt_vadv']['cmap']  = maps.jetWoGn()

defs['salt_hdiff']={}
defs['salt_hdiff']['var']   = 'salt_hdiff'
defs['salt_hdiff']['vmin']  =  -0.001
defs['salt_hdiff']['vmax']  =   0.001
defs['salt_hdiff']['label'] ='Salt horizontal diffusion [s$^ \mathrm{-1}$] '   
defs['salt_hdiff']['format']= '%.2g'
defs['salt_hdiff']['cmap']  = maps.jetWoGn()

defs['salt_vdiff']={}
defs['salt_vdiff']['var']   = 'salt_vdiff'
defs['salt_vdiff']['vmin']  =  -0.001
defs['salt_vdiff']['vmax']  =   0.001
defs['salt_vdiff']['label'] ='Salt vertical diffusion [s$^ \mathrm{-1}$] '   
defs['salt_vdiff']['format']= '%.2g'
defs['salt_vdiff']['cmap']  = maps.jetWoGn()


defs['salt_vflux']={}
defs['salt_vflux']['var']   = 'salt_vdiff'
defs['salt_vflux']['vmin']  =  -0.001
defs['salt_vflux']['vmax']  =   0.001
defs['salt_vflux']['label'] ='Salt vertical flux [s$^ \mathrm{-1}$] '   
defs['salt_vflux']['format']= '%.2g'
defs['salt_vflux']['cmap']  = maps.jetWoGn()


defs['salt_rate']={}
defs['salt_rate']['var']   = 'salt_rate'
defs['salt_rate']['vmin']  =  -0.001
defs['salt_rate']['vmax']  =   0.001
defs['salt_rate']['label'] ='Salt time rate change [s$^ \mathrm{-1}$] '   
defs['salt_rate']['format']= '%.2g'
defs['salt_rate']['cmap']  = maps.jetWoGn()

defs['dye']={}
defs['dye']['var']   = 'dye'
defs['dye']['vmin']  =  0.0
defs['dye']['vmax']  =  0.25
defs['dye']['label'] ='Dye concentration [kg m$^ \mathrm{-3}$] '   
defs['dye']['format']= '%.2g'
defs['dye']['cmap']  = maps.jetWoGn()

defs['dye_hadv']={}
defs['dye_hadv']['var']   = 'dye_hadv'
defs['dye_hadv']['vmin']  =  -0.0001
defs['dye_hadv']['vmax']  =   0.0001
defs['dye_hadv']['label'] ='Dye horizontal advection [kg m$^ \mathrm{-3}$ s$^ \mathrm{-1}$] '   
defs['dye_hadv']['format']= '%.2g'
defs['dye_hadv']['cmap']  = maps.jetWoGn()

defs['dye_vadv']={}
defs['dye_vadv']['var']   = 'dye_vadv'
defs['dye_vadv']['vmin']  =  -0.0001
defs['dye_vadv']['vmax']  =   0.0001
defs['dye_vadv']['label'] ='Dye vertical advection [kg m$^ \mathrm{-3}$ s$^ \mathrm{-1}$]  '   
defs['dye_vadv']['format']= '%.2g'
defs['dye_vadv']['cmap']  = maps.jetWoGn()

defs['dye_hdiff']={}
defs['dye_hdiff']['var']   = 'dye_hdiff'
defs['dye_hdiff']['vmin']  =  -0.0001
defs['dye_hdiff']['vmax']  =   0.0001
defs['dye_hdiff']['label'] ='Dye horizontal diffusion [kg m$^ \mathrm{-3}$ s$^ \mathrm{-1}$]  '   
defs['dye_hdiff']['format']= '%.2g'
defs['dye_hdiff']['cmap']  = maps.jetWoGn()

defs['dye_vdiff']={}
defs['dye_vdiff']['var']   = 'dye_vdiff'
defs['dye_vdiff']['vmin']  =  -0.0001
defs['dye_vdiff']['vmax']  =   0.0001
defs['dye_vdiff']['label'] ='Dye vertical diffusion [kg m$^ \mathrm{-3}$ s$^ \mathrm{-1}$]  '   
defs['dye_vdiff']['format']= '%.2g'
defs['dye_vdiff']['cmap']  = maps.jetWoGn()

defs['dye_vflux']={}
defs['dye_vflux']['var']   = 'dye_vdiff'
defs['dye_vflux']['vmin']  =  -0.0001
defs['dye_vflux']['vmax']  =   0.0001
defs['dye_vflux']['label'] ='Dye vertical diffusion [kg m$^ \mathrm{-3}$ s$^ \mathrm{-1}$]  '   
defs['dye_vflux']['format']= '%.2g'
defs['dye_vflux']['cmap']  = maps.jetWoGn()

defs['dye_rate']={}
defs['dye_rate']['var']   = 'dye_rate'
defs['dye_rate']['vmin']  =  -0.0001
defs['dye_rate']['vmax']  =   0.0001
defs['dye_rate']['label'] ='Dye time rate change [kg m$^ \mathrm{-3}$ s$^ \mathrm{-1}$]  '   
defs['dye_rate']['format']= '%.2g'
defs['dye_rate']['cmap']  = maps.jetWoGn()

defs['pmsl']={}
defs['pmsl']['var']   = 'pmsl'
defs['pmsl']['tvar']   = 'time'
defs['pmsl']['vmin']  =  0.
defs['pmsl']['vmax']  =  2.0e5
defs['pmsl']['label'] = 'pressure [Pa] '
defs['pmsl']['format']= '%.1g'
defs['pmsl']['cmap']  = maps.jetMinWi


## ADCIRC defs
defs['rad']={}
defs['rad']['var']   = 'radstress'
defs['rad']['tvar']   = 'time'
defs['rad']['vmin']  =  -1e-3
defs['rad']['vmax']  =   1e-3
defs['rad']['label'] = 'Wave force [m$^ \mathrm{-2}$ s$^ \mathrm{-2}$] '
defs['rad']['format']= '%.1g'
defs['rad']['cmap']  = maps.jetWoGn()



#plot extent
defs['lim']['xmin']  = -100.0 
defs['lim']['xmax']  =  400.0
defs['lim']['ymin']  =  1000.0
defs['lim']['ymax']  =  2000.0
defs['lim']['zmin']  = -6.0
defs['lim']['zmax']  =  1.0


