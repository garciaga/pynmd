from . import points_inside
from . import octant_plotting
from . import plot_settings
from . import colormaps
from . import vars_param
from . import taylor
#from . import plot_routines



# prepare cartopy BG env
import os
import string

path =  os.path.abspath(plot_settings.__file__).split('/')[:-1]
path = '/'.join(path)
os.environ["CARTOPY_USER_BACKGROUNDS"] = os.path.join(path,'BG/')
from . import general
