"""pvbatch --mesa sample_schlieren.py <flow.vtu> <out.npz>: density and d(rho)/dx (Gradient filter on the native mesh)
on a uniform (x, r) grid covering the schlieren field of view (nozzle exit -3.8 .. +35.2 mm) and the full radius of the
domain, for the line-of-sight projection. Also finds the nozzle-exit x (end of the lip wall at r 2.0-2.8 mm)."""
import os
import sys
import numpy as np
from paraview.simple import XMLUnstructuredGridReader, Gradient, ResampleToImage
from paraview import servermanager as sm
from vtk.util.numpy_support import vtk_to_numpy as v2n

vtu, out = sys.argv[1:3]
H = float(os.environ.get("H_MM", "0.03")) * 1e-3
rd = XMLUnstructuredGridReader(FileName=[vtu])
rd.UpdatePipeline()
b = rd.GetDataInformation().GetBounds()
print("domain bounds [mm]: x %.2f..%.2f  r %.2f..%.2f" % (b[0] * 1e3, b[1] * 1e3, b[2] * 1e3, b[3] * 1e3))


def sample(src, x0, x1, r0, r1, h):
    nx, nr = int(round((x1 - x0) / h)) + 1, int(round((r1 - r0) / h)) + 1
    rs = ResampleToImage(Input=src)
    rs.UseInputBounds = 0
    rs.SamplingBounds = [x0, x1, r0, r1, 0.0, 0.0]
    rs.SamplingDimensions = [nx, nr, 1]
    d = sm.Fetch(rs)
    pd = d.GetPointData()
    get = lambda n: v2n(pd.GetArray(n))
    ok = get("vtkValidPointMask").astype(bool).reshape(nr, nx)
    return np.linspace(x0, x1, nx), np.linspace(r0, r1, nr), ok, get, nr, nx


# nozzle exit: last solid x in the lip band
xs, rs_, ok, _, _, _ = sample(rd, 0.0, 4e-3, 2.0e-3, 2.8e-3, 0.005e-3)
xe = xs[np.where((~ok).any(axis=0))[0].max()]
print("nozzle exit x = %.4f mm" % (xe * 1e3))

g = Gradient(Input=rd)
g.ScalarArray = ["POINTS", "Density"]
g.ResultArrayName = "gradRho"
X0, X1 = xe - 3.8e-3, xe + 35.2e-3
R1 = b[3]
x, r, ok, get, nr, nx = sample(g, X0, X1, 0.0, R1, H)
rho = get("Density").astype(np.float32).reshape(nr, nx)
gr = get("gradRho").astype(np.float32).reshape(nr, nx, -1)[..., 0]
rho[~ok] = np.nan
gr[~ok] = np.nan
np.savez_compressed(out, x_mm=x * 1e3, r_mm=r * 1e3, rho=rho, drhodx=gr * 1e-3, valid=ok, xe_mm=xe * 1e3,
                    bounds_mm=np.array(b) * 1e3)
print("grid %d x %d, h %.3f mm, valid %d of %d; drho/dx in kg/m^3 per mm" % (nx, nr, H * 1e3, ok.sum(), ok.size))
print("saved", out)
