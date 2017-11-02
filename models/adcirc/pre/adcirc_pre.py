__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Portland State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"


import os




def fort14_update_bathy(data, path=None, file_name='fort.14'):
    """
    Write out bathymetry in data to ``fort.14`` formated file, by updating
    path/file_name accordingly

    :type data: dict(x,y,bathy) in vector format arranged based on node numbers
    :param data: python object to save the ``fort.14`` data to
    :type bathymetry: array or None
    :param bathymetry: if None then use the bathymetry in ``data``
    :type path: string or None
    :param path: path to the``fort.14`` fortmatted file
    :param string  file_name: file name

    """
    if path is None:
        path = os.getcwd()

    file_name = os.path.join(path, file_name)
    tmp = os.path.join(path, 'temp.14')

    # this currrently uses the present fort.14 as a template for formatting
    # purposes
    bathymetry = data['depth']
    x          = data['x']
    y          = data['y']
    nodes      = 1 + range(len(x))  
    # pylint: disable=C0103
    with open(file_name, 'r') as f, open(tmp, 'w') as fw:
        fw.write(f.readline())
        fw.write(f.readline())
        for i in range(len(x)):
            # pylint: disable=C0103
            fw.write('{:<7d} {:9.8E} {:9.8E} {:7.2f}\n'.format(nodes[i], x[i], y[i],
                                                               bathymetry[i]))
            f.readline()
        for line in f:
            fw.write(line)
    #os.rename(tmp, file_name)

