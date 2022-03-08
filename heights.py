import numpy as np
import math


def getHeights(M,facesID):
    faces_nodesID = M.faces.connectivities[facesID[:]] #return the nodes global ID
    adjacents_volumesID = M.faces.bridge_adjacencies(facesID[:],2,3)
    volume_center_coords = M.volumes.center(M.volumes())
    adjacents_volumes_coords = volume_center_coords[adjacents_volumesID[:,:]]

    face_nodes_coords = M.nodes.coords[faces_nodesID[:,0]] #one is enough
    normal_vector = M.faces.normal[facesID]
    normal_vector_modulus = np.linalg.norm(normal_vector,axis=1)
    d_coefficient = -(normal_vector*face_nodes_coords).sum(axis=1) # from the plane equation (a*x+b*y+c*z+d = 0)


    normal_vector_reshape = np.ones(adjacents_volumes_coords.shape)*normal_vector[:,np.newaxis,:]
    d_coefficient_reshape = np.ones(adjacents_volumesID.shape)*d_coefficient[:,np.newaxis]
    normal_vector_modulus_reshape = np.ones(adjacents_volumesID.shape)*normal_vector_modulus[:,np.newaxis]

    heights = abs((normal_vector_reshape*adjacents_volumes_coords).sum(axis=2)+d_coefficient_reshape)/normal_vector_modulus_reshape

    ''' OBS:
        normal_vector_reshape: reshaping the normal vector to the same shape as coords_vol
        d_coefficient_reshape: reshaping the d coefficient vector to the same shape as the hights
        normal_vector_modulus_reshape: reshaping the normal vector modulus term to the same shape as the hights
    '''
    return heights

def compute_heights(self):
    M = self.mesh
    internal_facesID = M.faces.internal[:]
    boundary_facesID = M.faces.boundary[:]
    internal_heights = getHeights(M,internal_facesID)
    boundary_heights = getHeights(M,boundary_facesID)

    return boundary_heights, internal_heights
#a*coords_vol[:,0]+b*coords_vol[:,1]+c*coords_vol[:,2]
