import unittest
import sys; sys.path.append("..")
import numpy as np
import plantbox as pb
from scipy.linalg import norm


def rootAge(l, r, k):  # root age at a certain length
    return -np.log(1 - l / k) * k / r


def rootLength(t, r, k):  # root length at a certain age
    return k * (1 - np.exp(-r * t / k))


def rootLateralLength(t, et, r, k):  # length of first order laterals (without second order laterals)
    i, l = 0, 0
    while et[i] < t:
        age = t - et[i]
        l += rootLength(age, r, k)
        i += 1
    return l


class TestRoot(unittest.TestCase):

    def root_example_rrp(self):
        """ an example used in the tests below, a main root with laterals """
        self.plant = pb.Organism()  # store organism (not owned by Organ, or OrganRandomParameter)
        p0 = pb.RootRandomParameter(self.plant)
        p0.name, p0.subType, p0.la, p0.lb, p0.lmax, p0.ln, p0.r, p0.dx = "taproot", 1, 10., 1., 100., 1., 1.5, 0.5
        p0.successor = [2]
        p0.successorP = [1.]
        p1 = pb.RootRandomParameter(self.plant)
        p1.name, p1.subType, p1.lmax, p1.r, p1.dx = "lateral", 2, 25., 2., 0.1
        self.p0, self.p1 = p0, p1  # needed at later point
        self.plant.setOrganRandomParameter(p0)  # the organism manages the type parameters and takes ownership
        self.plant.setOrganRandomParameter(p1)
        srp = pb.SeedRandomParameter(self.plant)
        self.plant.setOrganRandomParameter(srp)
        
        param0 = p0.realize()  # set up root by hand (without a root system)
        param0.la, param0.lb = 0, 0  # its important parent has zero length, otherwise creation times are messed up
        parentroot = pb.Root(1, param0, True, True, 0., 0., pb.Vector3d(0, 0, -1), 0, False, 0)  # takes ownership of param0
        parentroot.setOrganism(self.plant)
        parentroot.addNode(pb.Vector3d(0, 0, -3), 0)  # there is no nullptr in Python
        
        self.parentroot = parentroot  # store parent (not owned by child Organ)
        self.root = pb.Root(self.plant, p0.subType, pb.Vector3d(0, 0, -1), 0, self.parentroot , 0)
        self.root.setOrganism(self.plant)

    def test_dynamics(self):
        """ tests if nodes created in last time step are correct """  
        self.root_example_rrp()
        r = self.root
        r.simulate(.5, False)
        #self.assertEqual(r.hasMoved(), False, "dynamics: node is creaetd during first step")
        #r.simulate(1e-1, False)        
        #self.assertEqual(r.hasMoved(), True, "dynamics: node was expected to move, but did not")
        #non = r.getNumberOfNodes() 
        #r.simulate(2.4, False)
        #self.assertEqual(r.getOldNumberOfNodes(), non, "dynamics: wrong number of old nodes")
        #dx = r.getRootRandomParameter().dx
        #self.assertEqual(r.getNumberOfNodes() - non, round(2.4 * r.param().r / dx), "dynamics: unexpected number of new nodes")  # initially, close to linear growth




if __name__ == '__main__':
    unittest.main()
