from impress import FineScaleMeshMS as msh
from mpfa3d import mpfa3dScheme

M = msh('malha3D_Fluxo.msh', dim = 3)

tulio = mpfa3dScheme(M)
