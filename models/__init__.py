# Import all models
#
# Last Edit:
# 31 March 2015 - Gabriel Garcia Medina
from . import adcirc
from . import funwave
from . import nhwave
from . import roms
from . import swan
from . import ww3
from . import xbeach
from . import funwaveC
from . import parametric
from . import fvcom
from . import tools

try:
    from . import delft3d
except:
    print("Could not import Delft3D module")
    print("  check your okean installation and try again")
    print("  https://github.com/martalmeida/okean")
    
