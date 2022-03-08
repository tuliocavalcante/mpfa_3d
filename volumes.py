#volumes.py
import numpy as np

class Volumes:

    def __init__(self, pors, perm, mob, volume, nodes, faces, center):
        self.pors = pors
        self.perm = perm
        self.mob = mob
        self.nodes = nodes
        self.faces = faces
        self.center = center
        self.volume = volume

    def volume(self):
        v1 = self.nodes()
        v2 = self.nodes()
        v3 = self.nodes()
        volume = (1 / 6) * np.outer(v3,np.inner(v1, v2))
        return volume