from math import sqrt
#from Numeric import *
from numpy import *


k = 1000.0
dens = 997
d = 0.0011
pi = 3.14
vol = (4.0/3.0)*pi*(d/2.0)**3
m = vol*dens
v = array([2.12,0.102,0.5,0.0013,0.0076,0.00013,0.00084])
overlap = sqrt(m/k)*v
soverlap = overlap/d
print soverlap


