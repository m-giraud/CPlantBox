"""root system length over time"""
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import os
import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb

path = "../../../modelparameter/rootsystem/"
name = "Brassica_napus_a_Leitner_2010"

rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")
rs.initialize()

simtime = 200.  # days
dt = 1.
N = round(simtime / dt)  # steps

# Plot some scalar value over time
stype = "length"
v_, v1_, v2_, v3_ = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
for i in range(0, N):
    rs.simulate(dt)
    t = np.array(rs.getParameter("type"))
    v = np.array(rs.getParameter(stype))
    v_[i] = np.sum(v)
    v1_[i] = np.sum(v[t == 1])
    v2_[i] = np.sum(v[t == 2])
    v3_[i] = np.sum(v[t == 3])

t_ = np.linspace(dt, N * dt, N)
plt.plot(t_, v_, t_, v1_, t_, v2_, t_, v3_)
plt.xlabel("time (days)")
plt.ylabel(stype + " (cm)")
plt.legend(["total", "tap root", "lateral", "2. order lateral"])
plt.savefig("results/example_2d.png")
plt.show()

a = np.asarray([ t_, v_, v1_, v2_, v3_ ])  # !!!!
np.savetxt("temp.csv", np.transpose(a), delimiter=",")

# save final result 
ana = pb.SegmentAnalyser(rs)
ana.write("results/example_2d.vtp")

# Plot using paraview
os.system("cp results/example_2d.vtp temp.vtp")
os.system("paraview --script=../../pyscript/tubePlot.py")  # does not work, but works from termina???? works when run with pvpython

