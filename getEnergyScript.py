# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import numpy as np
import os
import sys

nGFS = 5000

ci = int(sys.argv[1])
Rhor = float(sys.argv[2])
Ohd = float(sys.argv[3])
Ohs = float(sys.argv[4])
Bond  = float(sys.argv[5])
We  = float(sys.argv[6])

name = "%d_getEnergy.dat" % ci

if os.path.exists(name):
    print("File %s found! New data will be appended to the file" % name)
for ti in range(nGFS):
    t = 0.010 * ti
    place = "intermediate/snapshot-%5.4f" % t
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
        exe = "./getEnergyAxi %s %s %s %s %s %s %s" % (place, name, Rhor, Ohd, Ohs, Bond, We)
        os.system(exe)
    print(("Done %d of %d" % (ti, nGFS)))
