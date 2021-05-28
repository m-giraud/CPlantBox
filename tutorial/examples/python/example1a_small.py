"""small example"""
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")

import plantbox as pb
import vtk_plot as vp

rs = pb.RootSystem()

# Open plant and root parameter from a file
path = "../../../modelparameter/rootsystem/"
name = "Brassica_oleracea_Vansteenkiste_2014"  # "Anagallis_femina_Leitner_2010" # 
rs.readParameters(path + name + ".xml")

# Initialize
rs.initialize()

# Simulate
rs.simulate(60, True)

# Export final result (as vtp)
rs.write("results/example_1a.vtp")
#
# # Plot, using vtk
vp.plot_roots(rs, "creationTime")

# Export results as segments (for animation)
ana = pb.SegmentAnalyser(rs)
ana.write("results/example_1a2.vtp")

# # Plot using paraview
# os.system("cp results/example_1a2.vtp temp.vtp")
# os.system("paraview --script=../../pyscript/tubePlot.py")  # does not work, but works from termina????l, works when run with pvpython

