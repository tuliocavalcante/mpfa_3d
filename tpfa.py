import numpy as np
import pdb

class tpfaScheme(object):
    def __init__(self, mesh):
        self.mesh = mesh
        self.node_flag, self.node_value, self.face_flag = self.calc_flagvalue(1)
        self.flag_value = {101:0, 102:1, 201:0, 202:-1, 203:1}
        self.permeabilities = {1: np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])}

    def __call__(self):
        self.run()
        pass

    def run(self):
        self.T, self.Q = self.assembly_tpfa()


    def assembly_tpfa(self):

        return T,Q

    def calc_flagvalue(self, caso):

        node_flag = 201*np.ones(len(self.mesh.nodes.all))
        face_flag = 201*np.ones(len(self.mesh.faces.all))
        coord = self.mesh.nodes.coords[:]
        node_value = np.zeros(len(self.mesh.nodes.all))

        for key in self.mesh.faces.flag:
            if key > 200:
                for i in self.mesh.faces.flag[key]:
                    face_flag[i] = key
                    noI = self.mesh.faces.bridge_adjacencies(self.mesh.faces.all, 0, 0)[i, 0]
                    noJ = self.mesh.faces.bridge_adjacencies(self.mesh.faces.all, 0, 0)[i, 1]
                    noK = self.mesh.faces.bridge_adjacencies(self.mesh.faces.all, 0, 0)[i, 2]
                    node_flag[noI] = key
                    node_flag[noJ] = key
                    node_flag[noK] = key

        for key in self.mesh.faces.flag:
            if key < 200:
                for i in self.mesh.faces.flag[key]:
                    face_flag[i] = key
                    noI = self.mesh.faces.bridge_adjacencies(self.mesh.faces.all, 0, 0)[i, 0]
                    noJ = self.mesh.faces.bridge_adjacencies(self.mesh.faces.all, 0, 0)[i, 1]
                    noK = self.mesh.faces.bridge_adjacencies(self.mesh.faces.all, 0, 0)[i, 2]
                    node_flag[noI] = key
                    node_flag[noJ] = key
                    node_flag[noK] = key

        for key in self.mesh.nodes.flag:
            for i in self.mesh.nodes.flag[key]:
                node_flag[i] = key

        for i in range(len(node_flag)):
            if node_flag[i] < 200:
                if caso==1:
                    node_value[i] = coord[i, 0]
                else:
                    node_value[i] = coord[i, 0]*coord[i, 0]

        return node_flag, node_value, face_flag

    # import pdb; pdb.set_trace()