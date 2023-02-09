#!/usr/bin/python
"""
Taking care of matplotlib versions for points_inside_poly function

adopted from:
 Module: plot_settings.py
 Author: Varun Hiremath <vh63@cornell.edu>
 Created: Thu,  2 Apr 2009 05:06:31 -0400

"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Oregon State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"

import pylab

# Symbols 
symbols = ['-','--','-.',':','.',',','o','^','v','<','>','s','+','x','D','d','1','2','3','4','h','H','p']
# Symbols + line
lps = [k+'-' for k in ['^','o','v','s','+','>','x','D','d','.','<','1','2','3','4','h','H','p']]

marker = ['o','^','s','D','v','x','>','d','.','x','<','+','1','2','3','4','h','H','p','o','v','s','>','D','d','.','x','<','+','1','2','3','4','h','H','p','o','v','s','>','D','d','.','x','<','+','1','2','3','4','h','H','p','o','v','s','>','D','d','.','x','<','+','1','2','3','4','h','H','p','o','v','s','>','D','d','.','x','<','+','1','2','3','4','h','H','p','o','v','s','>','D','d','.','x','<','+','1','2','3','4','h','H','p']
# Colors
colors = ['r','b','g','c','m','y','k','b','g','r','c','m','y','k','g','b','r','g','c','m','y','k','b','g','r','c','m','y','k','w','r','b','g','c','m','y','k','b','g','r','c','m','y','k','g','b','r','g','c','m','y','k','b','g','r','c','m','y','k','r','b','g','c','m','y','k','b','g','r','c','m','y','k','g','b','r','g','c','m','y','k','b','g','r','c','m','y','k']

cl = [    \
'red'    ,\
'blue'   ,\
'green'  ,\
'fuchsia',\
'grey'   ,\
'lime'   ,\
'maroon' ,\
'navy'   ,\
'olive'  ,\
'purple' ,\
'black'  ,\
'red'    ,\
'silver' ,\
'teal'   ,\
'aqua'   ,\
'yellow'  \
]

nice_colors = ["#348ABD", "#A60628","#7A68A6",  "#FF8C00", "#467821", "#CF4457", "#188487", "#E24A33", "#F4A460","#4B0082"]

ordinal=['a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)',
         'q)','r)','s)','t)','u)','v)','w)','x)','y)','z)','aa)','ab)','ac)','ad)','ae)',
         'af)','ag)','ah)','ai)','aj)','ak)','al)','am)','an)','ao)','ap)','aq)','ar)',
         'as)','at)','au)','av)','aw)','ax)','ay)','az)','ba)','bb)','bc)','bd)','be)',
         'bf)','bg)','bh)']

ordinal2=['I)','II)','III)','IV)','V)','VI)']

ordinal3=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u']

linewidth=[0.2,0.75,1.5,0.2,0.5,1,0.2,0.5,1,0.2,0.5,1,0.2,0.5,1,0.2,0.5,1,0.2,0.5,1,0.2,0.5,1,0.2,0.5,1]

mark_every=[5,9,13,8,12,16,5,9,13,8,12,16,5,9,13,8,12,16,5,9,13,8,12,16,5,9,13,8,12,16]

def get_figsize(fig_width_pt):
    inches_per_pt = 1.0/72.0                # Convert pt to inch
    golden_mean = (pylab.sqrt(5)-1.0)/2.0    # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]      # exact figsize
    return fig_size

# Publishable quality image settings for 2-column papers
params0 = {'backend': 'pdf',
          'axes.labelsize': 4,
          'text.fontsize': 4,
          'xtick.labelsize': 4,
          'ytick.labelsize': 4,
          #'legend.pad': 0.1,    # empty space around the legend box
          'legend.fontsize': 3,
          'lines.markersize': 2.5,
          'font.size': 4,
         # 'text.usetex': True,
          'figure.figsize': get_figsize(250)}

# Medium sized images
params1 = {'backend': 'pdf',
          'axes.labelsize': 8,
          'text.fontsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          #'legend.pad': 0.1,     # empty space around the legend box
          'legend.fontsize': 8,
          'lines.markersize': 4,
          'font.size': 8,
#          'text.usetex': True,
          'figure.figsize': get_figsize(500)}

# Large images (default)
params2 = {'backend': 'pdf',
          'axes.labelsize': 11,
          'text.fontsize': 11,
          'xtick.labelsize': 11,
          'ytick.labelsize': 11,
          #'legend.pad': 0.1,     # empty space around the legend box
          'legend.fontsize': 9,
           'lines.markersize': 3,
          'font.size': 11,
#          'text.usetex': True,
          'figure.figsize': get_figsize(800)}

def set_mode(mode):
    if mode == "publish":
        pylab.rcParams.update(params0)
    elif mode == "medium":
        pylab.rcParams.update(params1)
    else:
        pylab.rcParams.update(params2)

def set_figsize(fig_width_pt):
    pylab.rcParams['figure.figsize'] = get_figsize(fig_width_pt)

#pylab.rcParams.update(params2)


def find(l, s):
    for i in range(len(l)):
        if l[i].find(s)!=-1:
            return i
    return None # Or -1


def adjust_leg(leg):
    try:
        frame=leg.get_frame()
        frame.set_edgecolor('None')
        frame.set_facecolor('None')
    except:
        pass        

def smooth(x,beta):
    """ kaiser window smoothing """
    window_len=11
    # extending the data at beginning and at the end
    # to apply the window at the borders
    s = pylab.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    w = pylab.kaiser(window_len,beta)
    y = pylab.convolve(w/w.sum(),s,mode='valid')
    return y[5:len(y)-5]


def print_nc_dates(ref_grd,time_key):
    import netCDF4
    import netcdftime
    #time_key='ocean_time'
    nc=netCDF4.Dataset(ref_grd)
    ncv_obs=nc.variables
    utim_obs=netcdftime.utime(ncv_obs[time_key].units)
    sec_obs=ncv_obs[time_key][:]
    date_obs=utim_obs.num2date(sec_obs)
    print(date_obs)
    return date_obs



def mkd(t_vec,date_orig):
    """
    For plotting in decimal years, rather than datenum.

    t_vec : vector of matlab datenume (floats)
    date_orig : python date time of orig
    output: vector of number of days

    """
    t_orig = datetime2matlabdn(date_orig)
    return (t_vec - t_orig) / 365.25 + t_orig


def exp_func(x, a, b, c):
    return a*pylab.exp(-b*x)+c

def power_func(x, a, b,c):
    return a * x**b + c

def x2_func(x, a, b):
    return a * x**2 + b

# define a gaussian function to fit the data
def gaussian_func(x, a, b, c, d):
    val = a* pylab.exp(-(x - b)**2 / c**2) + d
    return val
    
    

# find the indices of the points in (x,y) closest to the points in (xi,yi)
def nearxy(x,y,xi,yi):
    ind=ones(len(xi),dtype=int)
    for i in arange(len(xi)):
        dist=sqrt((x-xi[i])**2+(y-yi[i])**2)
        ind[i]=dist.argmin()
    return ind
        
def replace_pattern_line(filename, pattern, line2replace):
    """
    replace the whole line if the pattern found
    
    """
    
    tmpfile = filename+'.tmp2'
    os.system(' cp  -f ' + filename + '  ' + tmpfile)
    tmp  = open(tmpfile,'r')
    fil  = open(filename,'w')
    for line in tmp:
        fil.write(line2replace if pattern in line else line)
        
    tmp.close()
    fil.close()
    os.system('rm ' + tmpfile  )  
    
    
def str_assing(inp_str,new_part,index_of_start):
    
    """
    replace part of a string
    
    """
    
    
    lenn = len(index_of_start)
    inp_str = inp_str[:index_of_start] + new_part + inp_str[index_of_start + lenn:]
    return inp_str    