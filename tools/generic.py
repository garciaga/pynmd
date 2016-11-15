"""
Generic python tools.
"""

from __future__ import division,print_function

def flatten_list_of_lists(x):
    """
    Code to flatten a list of lists
    
    SOURCE:
    -------
    http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
    """
    
    return [item for sublist in x for item in sublist]

