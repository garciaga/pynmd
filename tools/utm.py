import pyproj
import numpy as _np

_projections = {}


def zone(coordinates):
    if 56 <= coordinates[1] < 64 and 3 <= coordinates[0] < 12:
        return 32
    if 72 <= coordinates[1] < 84 and 0 <= coordinates[0] < 42:
        if coordinates[0] < 9:
            return 31
        elif coordinates[0] < 21:
            return 33
        elif coordinates[0] < 33:
            return 35
        return 37
    return int((coordinates[0] + 180) / 6) + 1


def letter(coordinates):
    return 'CDEFGHJKLMNPQRSTUVWXX'[int((coordinates[1] + 80) / 8)]


def project(coordinates,z=None,l=None):
    """
    Convert from spherical to UTM

    PARAMETERS:
    -----------
    coordiantes - Tuple of longitude and latitude (lon,lat)
    z           - (optional) the zone number
    l           - (optional) the zone letter

    RETURNS:
    --------
    z          - zone number
    l          - zone letter
    x          - easting
    y          - northing
    
    NOTES:
    ------
    1. If both z and l are passed as arguments to the function then the code 
       will use those for the transformation. Otherwise it will automatically
       detect the approapriate UTM zone. This is useful when dealind with 
       datasets across boundaries. 
    """
    if not z and not l:
        z = zone(coordinates)
        l = letter(coordinates)        
    if z not in _projections:
        _projections[z] = pyproj.Proj(proj='utm', zone=z, ellps='WGS84')
    x, y = _projections[z](coordinates[0], coordinates[1])
    if y < 0:
        y += 10000000
    return z, l, x, y


def unproject(z, l, x, y):
    """
    Convert from UTM to spherical

    PARAMETERS:
    -----------
      z - zone Number
      l - Zone letter
      x - easting [m]
      y - northing [m]

    RETURNS:
    --------
      lon - longitude
      lat - latitude
    """
    if z not in _projections:
        _projections[z] = pyproj.Proj(proj='utm', zone=z, ellps='WGS84')
    if l < 'N':
        y -= 10000000
    lng, lat = _projections[z](x, y, inverse=True)
    return (lng, lat)

def haversine(lat1,lon1,lat2,lon2):
    """
    Compute distance between two points using the haversine formula

    RETURNS:
    --------
    d : distance between lat1,lon1 and lat2,lon2 in meters

    """
    sin_dlon = _np.sin(_np.pi / 180 * (lon2 - lon1) / 2)**2
    sin_dlat = _np.sin(_np.pi / 180 * (lat2 - lat1) / 2)**2
    x = sin_dlat + (_np.cos(_np.pi / 180 * lat1) * _np.cos(_np.pi / 180 * lat2) *
                    sin_dlon)
    d = 2 * 6378137 * _np.arcsin((x)**0.5)
    return d