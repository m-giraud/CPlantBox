""" map root segments to a soil grid """
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import numpy as np
import plantbox as pb
import vtk_plot as vp

""" parameters """
sim_time = 14  # [day]
rs_age = 1  # initial age
dt = 0.05  # [days] Time step must be very small
periodic = False

""" root system """
rs = pb.MappedRootSystem()
path = "../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml") 

""" soil """
min_ = np.array([-5, -5, -20])
max_ = np.array([5, 5, 0.])
res_ = np.array([1, 5, 10])
if not periodic:
    sdf = pb.SDF_PlantBox(0.99 * (max_[0] - min_[0]), 0.99 * (max_[1] - min_[1]), 0.99 * (max_[2] - min_[2]))
    rs.setGeometry(sdf)
rs.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments

rs.getRootSystemParameter().seedPos = pb.Vector3d(0., 0., -0.1)
rs.setSeed(100)
rs.initialize()
rs.simulate(rs_age, False)
N = round(sim_time / dt)

ana = pb.SegmentAnalyser(rs.mappedSegments())
anim = vp.AnimateRoots(ana)
anim.min = min_
anim.max = max_
anim.res = res_
anim.avi_name = "avi/example_7a_"
anim.start()

for i in range(0, N):

    rs.simulate(dt, False)

    """ add segment indices """
    segs = rs.segments
    x = np.zeros(len(segs))
    for i, s in enumerate(segs):
        try:
            x[i] = rs.seg2cell[i]
        except:  # in case the segment is not within the domain
            x[i] = -10

    ana = pb.SegmentAnalyser(rs.mappedSegments())
    ana.addData("linear_index", x)

    anim.rootsystem = ana
    anim.root_name = "linear_index"
    anim.update()
