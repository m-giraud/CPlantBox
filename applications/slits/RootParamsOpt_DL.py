"""Spring barley root parameter optimization with std"""
import sys
sys.path.append("../..")
sys.path.append("../../src/python_modules")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vtk_plot as vp
import plantbox as pb
import pickle as pkl
import math
import numpy.matlib
from scipy.optimize import differential_evolution
from pyevtk.hl import gridToVTK  # pip3 install pyevtk
import time

path = "../../../modelparameter/rootsystem/"
name = "spring_barley_Fernand"

t = time.process_time()

simtime = 107
N = 24  # number of rows, 8 pro [m]
M = 2  # number of plants in a row 
distp = 3.  # distance between the root systems along row[cm]
distr = 12.5  # distance between the rows[cm]
distTr = N * distr  # total row spacing
distTp = M * distp  # total distance between plants 
top = 0.  # vertical top position (cm) (normally = 0)
bot = -100.  # vertical bot position (cm) (e.g. = -100 cm)
left = -147.5  # left along y-axis (cm)
right = 147.5  # right along y-axis (cm)
n = 20  # number of layers, each with a height of (top-bot)/n
m = 60  # number of horizontal grid elements (each with length of (right-left)/m)
m2 = int(M * distp)  # resolution of 1 cm  
exact = True  # calculates the intersection with the layer boundaries (true), only based on segment midpoints (false)

soilVolume = 5 * 5 * 5

times = [73, 72, 71, 46, 45]

# Measured RLD
with open('RLD_DL.pkl', 'rb') as f:
        measured_RLD = pkl.load(f)


def err(fitparams):
        lmax0 		 = fitparams[0]
        theta0   	 = fitparams[1]
        r0 		 = fitparams[2]
        ln0 		 = fitparams[3]
        tropismN0 	 = fitparams[4]
        lb1 		 = fitparams[5]
        la1 		 = fitparams[6]
        lmax1 		 = fitparams[7]
        maxB0 		 = fitparams[8]

        # Define bulk denistitis in slit
        # print("EquidistantGrid3D: ", right - left, distp * M, top - bot, "cm, res:", m, m2, n)
        scale_elongation = pb.EquidistantGrid3D(right - left, distp * M, top - bot, m, m2, n)  # grid is centered in x, y 
        bulk_density = np.ones((m, m2, n))

        # all bulk values
        zi0 = n - 1 - int(30 / ((top - bot) / n))  # gridz is from [-depth - 0], i.e. surface is at end of array
        zi1 = n - 1 - int(50 / ((top - bot) / n))
        # print("z indices", zi0, zi1)
        bulk_density[:,:, zi0:] = 1.7  # 0 - 30 cm 
        bulk_density[:,:, zi1:zi0] = 1.786  # 30 - 50 cm
        bulk_density[:,:, 0:zi1] = 1.786 + 0.1  # 50 - 100 cm (? + 0.1 for vizualisation)

        # slit only
        xi1 = int(135 / ((right - left) / m))  # -15 cm
        xi2 = int(165 / ((right - left) / m))  # 15 cm
        # print("x indices", xi1, xi2)
        bulk_density[xi1:xi2,:, zi0:] = 1.26  # -150 to -15 cm 
        bulk_density[xi1:xi2,:, zi1:zi0] = 1.4  # -15 to 15 cm
        bulk_density[xi1:xi2,:, 0:zi1] = 1.575 + 0.2  # 15 - 150 cm (?  +0.2 for vizualisation)

        scales = 39.553 * np.exp(-1.851 * bulk_density) / 1.68  #  equation TODO, see Mondrage et al.
        scale_elongation.data = scales.flatten('F')  # set proportionality factors
        
        allRS = []
        for i in range(0, N):
                for j in range(0, M):
                        rs = pb.RootSystem()
                        rs.readParameters(path + name + ".xml")
                        for p in rs.getRootRandomParameter():  # set scale elongation function for all root types
                                p.f_se = scale_elongation
                        p0 = rs.getRootRandomParameter(1)  # tap and basal root type
                        p1 = rs.getRootRandomParameter(2)
                        srp = rs.getRootSystemParameter()

                        p0.lmax 	 = lmax0
                        p0.theta 	 = theta0
                        p0.r 		 = r0
                        p0.ln 	 = ln0
                        p0.tropismN 	 = tropismN0
                        p1.lb   	 = lb1
                        p1.la 	 = la1
                        p1.lmax 	 = lmax1
                        srp.maxB 	 = maxB0

                        rs.setSeed(1)    
                        rs.getRootSystemParameter().seedPos = pb.Vector3d(left + distr * (i + 0.5), -distp / 2 * M + distp * (j + 0.5), -3.)  # cm
                        rs.initialize(False)
                        allRS.append(rs)

        # # Simulate
        # time = 0
        # dt = 1                                                                         # day
        # while time < simtime:                                                          # for future coupling with dynamic water movement 
        #         # print("day", time)     

        # update scales (e.g. from water content, soil_strength)
        scales = 39.553 * np.exp(-1.851 * bulk_density) / 1.68  # equation TODO, see Mondrage et al.
        scale_elongation.data = scales.flatten('F')
        
        for rs in allRS:
                rs.simulate(simtime)
        # time += dt

        # Export results as single vtp files (as polylines)
        ana = pb.SegmentAnalyser()  # see example 3b
        for z, rs in enumerate(allRS):
                # vtpname = "results/plantsb" + str(z) + ".vtp"
                # rs.write(vtpname)
                ana.addSegments(rs)  # collect all

        # Write all into single file (segments)
        # ana.write("results/plantsb_allDL.vtp")

        # Set periodic domain
        ana.mapPeriodic(distTr, distTp)     
        ana.pack()                      
        # ana.write("results/plantsb_periodicDL.vtp")

        rl_ = []
        for j in range(len(times)):
                # print("creating rld for time", times[j])
                ana.filter("creationTime", 0, times[j])
                rld = ana.distribution2("length", top, bot, left, right, n, m, False)
                rl_.append(rld)
                
        rl_.reverse()
        rld_ = np.array(rl_) / soilVolume

        df = rld_.reshape(-1, rld_.shape[1])
                
        np.savetxt("sim_dataoptDL.txt", np.round(df, decimals=4), fmt='%.4f', delimiter=',')

        # reduce dimension to 100x20
        df2 = pd.DataFrame(df)
        lst = [i for i in range(20, 60)] + [i for i in range(80, 120)] + [i for i in range(140, 180)] \
              +[i for i in range(200, 240)] + [i for i in range(260, 300)]
        rld_simDL = df2.drop(lst)
                
        # print("dimension of simDL:", rld_simDL.shape)

        err = np.sqrt(sum((rld_simDL - measured_RLD) ** 2) / len(measured_RLD))
        return err
                

bounds = ([100, 200], [1.5, 2], [1.5, 3], [0, 2], [1, 5], [0, 2], [0, 3], [0, 10], [5, 20])  # [lmax0,theta0,r0,ln0,tropismN0,lb1,la1,lmax1,maxB0]

result = differential_evolution(err, bounds, strategy='best2bin', maxiter=1000, seed=3, popsize=20, mutation=0.5) 
print(result.x)

elapsed_time = (time.process_time() - t) / 3600

print("Time for program to execute in hours: ", elapsed_time)

# lmax0  = 
# theta0 = 
# r0     = 
# ln0    = 
# tropismN0 =                             
# lb1    = 
# la1    =
# lmax1  = 
# maxB0  = 
# time to excute sim = 

