from impress import FineScaleMeshMS as msh
from mpfa3dcorea import mpfa3dSchemecorea

M = msh('malha3D_Fluxo.msh', dim = 3)

tulio = mpfa3dSchemecorea(M)
