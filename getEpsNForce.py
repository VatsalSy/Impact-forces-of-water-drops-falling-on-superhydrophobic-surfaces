# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import numpy as np
import os
import sys

nGFS = 27510

ci = int(sys.argv[1])
Ohd = float(sys.argv[2])
We  = float(sys.argv[3])

folder = 'bview'  # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)

name = "%4.4d_EpsForce.dat" % ci

if os.path.exists(name):
    print("File %s found! New data will be appended to the file" % name)
for ti in range(nGFS):
    t = 0.01 * ti
    place = "intermediate/snapshot-%5.4f" % t
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
        exe = "./getEpsForce %s %s %s %s" % (place, name, Ohd, We)
        os.system(exe)
    print(("Done %d of %d" % (ti, nGFS)))
