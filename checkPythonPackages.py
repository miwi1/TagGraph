print ("")

try:
    import sqlalchemy
    print ('sqlalchemy exist')
except ImportError:
    print ('please install sqlalchemy: pip install sqlalchemy')

try:
    import pympler
    print ('pympler exist')
except ImportError:
    print ('please install pympler: pip install pympler')

try:
    import numpy
    print ('numpy exist')
except ImportError:
    print ('please install numpy: yum install numpy')

try:
    import networkx
    print ('networkx exist')
except ImportError:
    print ('please install networkx: pip install networkx')
		
try:
    import pymzml
    print ('pymzml exist')
except ImportError:
    print ('please install pymzml: pip install pymzml')

try:
    import MySQLdb
    print ('MySQLdb exist')
except ImportError:
    print ('please install MySQLdb: yum install MySQL-python')

try:
    import pyteomics
    print ('pyteomics exist')
except ImportError:
    print ('please install pyteomics: pip install pyteomics')

print ("")
