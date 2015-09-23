import matplotlib


"""
Taking care of matplotlib versions for points_inside_poly function
"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2013, Oregon State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"


if (matplotlib.__version__ < '1.2'):
    from matplotlib.nxutils import points_inside_poly
else:
    from matplotlib.path import Path

def inside_poly(data,vertices):
  if(matplotlib.__version__ < '1.2'):
      mask = points_inside_poly(data, vertices)
  else:
      mask = Path(vertices).contains_points(data)
  return mask


#use like this
#from points_inside import inside_poly as points_inside_poly
#
#
#
#
