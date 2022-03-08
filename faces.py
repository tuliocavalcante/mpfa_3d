#faces.py
import numpy as np

class Faces:

    def __init__(self, perm, mob, area, nodes, erel):
        self.nodes = nodes
        self.erel = erel
        self.area = area
        self.perm = perm
        self.mob = mob

    def area(self):
        v1 = self.nodes()
        v2 = self.nodes()
        volume = (1 / 2) * np.outer(v1, v2)
        return volume

    def perm(self):
        p1 = self.erel()
        p2 = self.erel()
        return perm