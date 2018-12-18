# Taken from http://www.mobify.com/blog/sqlalchemy-memory-magic/
from pympler import summary, muppy
#import psutil

class StringFolder(object):
    """
    Class that will fold strings. See 'fold_string'.
    This object may be safely deleted or go out of scope when
    strings have been folded.
    """
    def __init__(self):
        self.unicode_map = {}
            
    def fold_string(self, s):
        """
        Given a string (or unicode) parameter s, return a string object
        that has the same value as s (and may be s). For all objects
        with a given value, the same object will be returned. For unicode
        objects that can be coerced to a string with the same value, a
        string object will be returned.
        If s is not a string or unicode object, it is returned unchanged.
        :param s: a string or unicode object.
        :return: a string or unicode object.
        """
        # If s is not a string or unicode object, return it unchanged
        if not isinstance(s, basestring):
            return s
                                                    
        # If s is already a string, then str() has no effect.
        # If s is Unicode, try and encode as a string and use intern.
        # If s is Unicode and can't be encoded as a string, this try
        # will raise a UnicodeEncodeError.
        try:
            return intern(str(s))
        except UnicodeEncodeError:
            # Fall through and handle s as Unicode
            pass
        
        # Look up the unicode value in the map and return
        # the object from the map. If there is no matching entry,
        # store this unicode object in the map and return it.
        t = self.unicode_map.get(s, None)
        if t is None:
            # Put s in the map
            t = self.unicode_map[s] = s
        return t

folder = StringFolder()

def string_folding_wrapper(results):
    """
    This generator yields rows from the results as tuples,
    with all string values folded.
    """
    # Get the list of keys so that we build tuples with all
    # the values in key order.
    keys = results.keys()
    for row in results:
        yield tuple( folder.fold_string(row[key]) for key in keys )


def fold_list(item):
    return tuple(folder.fold_string(element) for element in item)

'''
# This doesn't belong here, but it's going here anyways
def get_virtual_memory_usage_kb():
    """
    The process's current virtual memory size in Kb, as a float.
    
    """
    return float(psutil.Process().memory_info_ex().vms) / 1024.0

def memory_usage(where):
    """
    Print out a basic summary of memory usage.
    
    """
    mem_summary = summary.summarize(muppy.get_objects())
    print "Memory summary:", where
    summary.print_(mem_summary, limit=2)
    print "VM: %.2fMb" % (get_virtual_memory_usage_kb() / 1024.0)
'''
