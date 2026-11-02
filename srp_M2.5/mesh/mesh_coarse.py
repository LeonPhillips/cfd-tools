"""Coarse axisymmetric RANS mesh for srp_M2.5 (supersonic retropropulsion, M = 2.5).

Confirmed assumptions (2026-10-26): units metres; freestream flows in -x from edge 16, jet
exits in +x; axis y = 0; air, ideal gas, Sutherland; adiabatic no-slip walls; wall-resolved
y+ ~ 1, no wall functions. Edge -> marker map is in MARKERS below (checked against the user's
labelled screenshots).
v2 (2026-10-27): first cell 1.8e-7 m on all walls, sized for the nozzle throat at
p0j = 8.015 bar (user approved); v1 used 3.875e-7 m from the freestream at p0inf = 6 bar.
v3 (2026-10-27): 3 mm cells from the end of the near-wake box (x = -300 mm) to the outlet, full
height, and the 1.5 mm wake box ramps out at 0.05 per mm instead of 0.2 (user approved). In v2 the
cells jumped 1.5 -> 10 mm over 42.5 mm there, and the implicit solver became unstable at CFL ~5.8.
v4 (2026-10-27): new geometry su2_shortDiff_wFillet_v4.STEP (user): the duct ends at a vertical outlet at
x = -308.1 mm (r 0..55 mm), the base corner fillet is 0.5 mm (was 0.1 mm) and the diverging slip wall rises
r 45 -> 55 mm over 250 mm. Curve numbers changed (marker map approved by the user, figure
edge_map_v4.png). The 1.5 mm near-wake box now reaches the outlet, so there is no size transition in the wake.
meshMedium (2026-10-28): optional MED_BOX env var adds a finer box; unset, the mesh is exactly v4 (meshCoarse).
Usage: python mesh_coarse.py <step file> <output .su2>
"""
import math
import os
import sys

import gmsh

STEP, OUT = sys.argv[1], sys.argv[2]

# ---- flow conditions used only to size the wall layers ---------------------------------
GAMMA, R_GAS = 1.4, 287.058
MACH, P0, T0 = 2.5, 3.247e5, 288.15   # freestream stagnation (case SA_p0j8.015_p0inf3.247)
P0J = 8.015e5                          # jet stagnation pressure
X_REF = 0.018            # running length for the skin-friction estimate (body diameter)
YPLUS = 1.0


def sutherland(t):
    return 1.716e-5 * (t / 273.15) ** 1.5 * (273.15 + 110.4) / (t + 110.4)


def first_cell(rho, u, mu, x=X_REF, yplus=YPLUS):
    re_x = rho * u * x / mu
    cf = 0.026 / re_x ** (1 / 7)            # turbulent flat plate (White)
    u_tau = math.sqrt(0.5 * cf) * u
    return yplus * mu / (rho * u_tau), re_x


ratio = 1 + 0.5 * (GAMMA - 1) * MACH ** 2
t_inf, p_inf = T0 / ratio, P0 / ratio ** (GAMMA / (GAMMA - 1))
rho_inf = p_inf / (R_GAS * t_inf)
u_inf = MACH * math.sqrt(GAMMA * R_GAS * t_inf)
mu_inf = sutherland(t_inf)
h_body, re_x = first_cell(rho_inf, u_inf, mu_inf)
print(f"freestream: T={t_inf:.2f} K  p={p_inf:.1f} Pa  rho={rho_inf:.4f} kg/m3  "
      f"U={u_inf:.1f} m/s  mu={mu_inf:.3e}  Re/m={rho_inf*u_inf/mu_inf:.3e}")
print(f"body wall, y+={YPLUS} (x_ref={X_REF} m, Re_x={re_x:.3e}): {h_body:.3e} m")

# Nozzle throat (sonic) at p0j: the most demanding wall, so it sets the first cell everywhere
t_s, rho_s = T0 / 1.2, P0J / (R_GAS * T0) * 1.2 ** (-1 / (GAMMA - 1))
u_s = math.sqrt(GAMMA * R_GAS * t_s)
h_noz, _ = first_cell(rho_s, u_s, sutherland(t_s), x=0.004)
print(f"nozzle throat at p0j={P0J/1e5:.3f} bar, y+={YPLUS}: {h_noz:.3e} m")
H1 = float(os.environ.get("H1", "1.8e-7"))
print(f"first cell height used: {H1:.3e} m")

# ---- mesh parameters (metres) ----------------------------------------------------------
BL_RATIO = 1.2
# v1 (H1 3.875e-7): 16 layers = 0.034 mm; 17 meshed, 18 failed ("2 intersections").
# v2 (H1 1.8e-7): 20 layers = 0.0336 mm, the same total thickness as v1.
BL_LAYERS = int(os.environ.get("BL_LAYERS", "20"))
BL_THICK = H1 * (BL_RATIO ** BL_LAYERS - 1) / (BL_RATIO - 1)
S_FILLET = 1.0e-5        # 0.01 mm on the 0.03 mm fillets
S_WALL = 2.0e-5          # 0.02 mm next to all no-slip walls (top of boundary layer)
S_NOZZLE = 5.0e-5        # nozzle interior
S_REFINE = S_NOZZLE      # plume refinement box x 0..30 mm, r 0..9 mm (user, 2026-10-26)
REFINE_GROWTH = 0.1      # growth away from the plume box (user, 2026-10-26)
S_JET = 6.0e-4           # max size in the jet / bow-shock region (was 0.4 mm; user set 0.6 mm)
S_BODY = 6.0e-4          # around the body side
S_WAKE = 1.5e-3          # near wake
WAKE_GROWTH = 0.05       # ramp out of the near-wake box (v3; was GROWTH = 0.2)
S_FAR = 1.0e-2           # rest of the domain

# v4 geometry (su2_shortDiff_wFillet_v4.STEP)
MARKERS = {
    "axis": [4, 9],
    "inlet_freestream": [5],
    "inlet_jet": [3],
    "wall_noslip": [1, 2] + list(range(10, 20)),
    "wall_slip": [6, 7],
    "outlet": [8],
}
FILLETS = [1, 11, 13, 16, 18]  # 0.1 mm fillets, plus the 0.5 mm base-corner fillet (curve 11)

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setString("Geometry.OCCTargetUnit", "M")
gmsh.model.add(os.path.basename(OUT).replace(".su2", ""))
gmsh.model.occ.importShapes(STEP)
gmsh.model.occ.synchronize()

curves = sorted(t for _, t in gmsh.model.getEntities(1))
assert curves == list(range(1, 20)), f"unexpected curve tags {curves}"
for name, tags in MARKERS.items():
    gmsh.model.addPhysicalGroup(1, tags, name=name)
gmsh.model.addPhysicalGroup(2, [t for _, t in gmsh.model.getEntities(2)], name="fluid")

f = gmsh.model.mesh.field
walls = MARKERS["wall_noslip"]

USE_BL = os.environ.get("USE_BL", "1") == "1"
ALGO = int(os.environ.get("MESH_ALGO", "8"))
print(f"variant: USE_BL={USE_BL} MESH_ALGO={ALGO}", flush=True)

# Wall-layer fields, one per group (BL_GROUPS env var: groups separated by ";", each a comma list of
# curves or ranges, e.g. "1-2,10-19"). A single group is used: splitting did not avoid the 18-layer
# intersection; only a thinner layer did.
def parse_groups(spec):
    groups = []
    for g in spec.split(";"):
        group = []
        for item in g.split(","):
            a, _, b = item.partition("-")
            group += list(range(int(a), int(b or a) + 1))
        groups.append(group)
    return groups


BL_GROUPS = parse_groups(os.environ.get("BL_GROUPS", "1-2,10-19"))


def end_points(c1, c2):
    """Geometric point shared by curves c1 and c2."""
    b = lambda c: {t for _, t in gmsh.model.getBoundary([(1, c)], oriented=False)}
    return sorted(b(c1) & b(c2))


# Points where the wall layers meet a non-wall boundary at 90 degrees: base wall (10) / axis (9)
# and plenum wall (2) / jet inlet (3). Listing them in PointsList makes the layers end cleanly
# against that boundary; without it Gmsh stopped short and filled the gap with sliver triangles.
BL_END_POINTS = end_points(9, 10) + end_points(2, 3)
if USE_BL:
    print(f"wall-layer end points: {BL_END_POINTS}", flush=True)
    print(f"wall-layer groups: {BL_GROUPS}", flush=True)
    for group in BL_GROUPS:
        bl = f.add("BoundaryLayer")
        f.setNumbers(bl, "CurvesList", group)
        f.setNumber(bl, "Size", H1)
        f.setNumber(bl, "Ratio", BL_RATIO)
        f.setNumber(bl, "Thickness", BL_THICK)
        f.setNumber(bl, "Quads", 1)
        if os.environ.get("BL_USE_END_POINTS", "1") == "1":
            f.setNumbers(bl, "PointsList", BL_END_POINTS)
        if "BL_ANISOMAX" in os.environ:
            f.setNumber(bl, "AnisoMax", float(os.environ["BL_ANISOMAX"]))
        f.setAsBoundaryLayer(bl)


def threshold(curve_list, s_min, s_max, d_min, d_max, local=False):
    d = f.add("Distance")
    f.setNumbers(d, "CurvesList", curve_list)
    f.setNumber(d, "Sampling", 400)
    t = f.add("Threshold")
    f.setNumber(t, "InField", d)
    f.setNumber(t, "SizeMin", s_min)
    f.setNumber(t, "SizeMax", s_max)
    f.setNumber(t, "DistMin", d_min)
    f.setNumber(t, "DistMax", d_max)
    # local=True: no effect beyond DistMax (otherwise SizeMax applies everywhere and
    # wins in the Min field)
    f.setNumber(t, "StopAtDistMax", 1 if local else 0)
    return t


# Cell size may grow by at most GROWTH per unit distance away from walls and boxes
# (0.2 -> neighbouring cells differ by ~20 %). Steeper ramps left 0.02 mm cells next to
# >1 mm cells and produced collapsed slivers on the axis behind the base.
GROWTH = 0.2


def box(x0, x1, y0, y1, size, v_out=S_FAR, growth=GROWTH):
    b = f.add("Box")
    for k, v in dict(XMin=x0, XMax=x1, YMin=y0, YMax=y1, ZMin=-1, ZMax=1,
                     VIn=size, VOut=v_out, Thickness=(v_out - size) / growth).items():
        f.setNumber(b, k, v)
    return b


fields = [
    threshold(FILLETS, S_FILLET, S_WALL, 2e-5, 2e-4, local=True),
    threshold(walls, S_WALL, S_FAR, BL_THICK, BL_THICK + (S_FAR - S_WALL) / GROWTH),
    box(-8.2e-3, 2.0e-3, 0.0, 3.0e-3, S_NOZZLE),               # nozzle + plenum (VOut applies everywhere outside, so keep S_FAR)
    box(0.0, 30e-3, 0.0, 9e-3, S_REFINE, growth=REFINE_GROWTH),  # plume refinement
    box(-10e-3, 75e-3, 0.0, 46e-3, S_JET),                     # jet / bow-shock region (cap)
    box(-60e-3, -10e-3, 0.0, 20e-3, S_BODY),                   # body side
    box(-310e-3, -10e-3, 0.0, 60e-3, S_WAKE, growth=WAKE_GROWTH),  # near wake, to the outlet (x -308.1 mm)
]
# Optional extra refinement boxes (MED_BOX="x0,x1,r0,r1,size" in mm; several boxes separated by ";"). meshMedium
# (user, 2026-10-28): x 0.5..20 mm, r 0..9.5 mm at 0.025 mm, around the jet, Mach disk and bow shock.
# meshFineNozzle (user, 2026-10-28): that box at 0.0177 mm plus the nozzle + plenum x -8.2..0.5 mm, r 0..3 mm at 0.0177 mm.
MED_BOX = os.environ.get("MED_BOX")
for spec in filter(None, (MED_BOX or "").split(";")):
    bx0, bx1, by0, by1, bsize = (float(v) * 1e-3 for v in spec.split(","))
    fields.append(box(bx0, bx1, by0, by1, bsize, growth=REFINE_GROWTH))
    print(f"extra refinement box: x {bx0*1e3:g}..{bx1*1e3:g} mm, r {by0*1e3:g}..{by1*1e3:g} mm, size {bsize*1e3:g} mm",
          flush=True)
fmin = f.add("Min")
f.setNumbers(fmin, "FieldsList", fields)
f.setAsBackgroundMesh(fmin)

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.Algorithm", ALGO)            # 8 frontal-Delaunay for quads, 6 frontal-Delaunay
# 0 = simple, 1 = blossom (crashes with the boundary-layer field in Gmsh 4.15.2), -1 = off
RECOMB = int(os.environ.get("RECOMB_ALGO", "0"))
print(f"recombination: {RECOMB}", flush=True)
gmsh.option.setNumber("Mesh.RecombineAll", 0 if RECOMB < 0 else 1)
if RECOMB >= 0:
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", RECOMB)
gmsh.option.setNumber("Mesh.Smoothing", 5)

gmsh.model.mesh.generate(2)

types, tags, _ = gmsh.model.mesh.getElements(2)
counts = {gmsh.model.mesh.getElementProperties(t)[0]: len(e) for t, e in zip(types, tags)}
print("2D elements:", counts, "total", sum(counts.values()))
print("nodes:", len(gmsh.model.mesh.getNodes()[0]))
all2d = [e for tl in tags for e in tl]
q = gmsh.model.mesh.getElementQualities(all2d, "minSICN")
print(f"minSICN: min={min(q):.3f}  <0.1: {sum(x < 0.1 for x in q)}  <0: {sum(x <= 0 for x in q)}")
print(f"BL: first cell {H1:.3e} m, {BL_LAYERS} layers, ratio {BL_RATIO}, thickness {BL_THICK:.3e} m")

gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(OUT)
gmsh.write(OUT.replace(".su2", ".msh"))


def write_vtk_files(base):
    """ParaView files (legacy ASCII VTK, coordinates in metres):
    <base>.vtk            2D cells with cell data: quality (Gmsh minSICN), size (sqrt area), cell type
    <base>_boundaries.vtk boundary edges with cell data marker_id; names in the header comment
    """
    ntags, coords, _ = gmsh.model.mesh.getNodes()
    index = {t: i for i, t in enumerate(ntags)}
    xyz = coords.reshape(-1, 3)
    points = "\n".join(f"{x:.9e} {y:.9e} 0" for x, y, _ in xyz)
    vtk_type = {3: 5, 4: 9, 2: 3}    # nodes per element -> VTK triangle, quad, line

    def cells_of(dim, tag=-1):
        out = []
        for etype, etags, enodes in zip(*gmsh.model.mesh.getElements(dim, tag)):
            n = gmsh.model.mesh.getElementProperties(etype)[3]
            out += [(t, [index[k] for k in e]) for t, e in zip(etags, enodes.reshape(-1, n))]
        return out

    def write(path, cells, title, data):
        size = sum(len(c) + 1 for c in cells)
        with open(path, "w") as fh:
            fh.write(f"# vtk DataFile Version 3.0\n{title}\nASCII\nDATASET UNSTRUCTURED_GRID\n")
            fh.write(f"POINTS {len(xyz)} double\n{points}\n")
            fh.write(f"CELLS {len(cells)} {size}\n")
            fh.write("\n".join(f"{len(c)} " + " ".join(map(str, c)) for c in cells) + "\n")
            fh.write(f"CELL_TYPES {len(cells)}\n" + "\n".join(str(vtk_type[len(c)]) for c in cells) + "\n")
            fh.write(f"CELL_DATA {len(cells)}\n")
            for name, kind, values in data:
                fh.write(f"SCALARS {name} {kind} 1\nLOOKUP_TABLE default\n")
                fh.write("\n".join(str(v) for v in values) + "\n")
        print("wrote", path)

    surf = cells_of(2)
    quality = gmsh.model.mesh.getElementQualities([t for t, _ in surf], "minSICN")
    xy = xyz[:, :2]
    size = [math.sqrt(0.5 * abs(sum(xy[c[i]][0] * xy[c[(i + 1) % len(c)]][1] -
                                    xy[c[(i + 1) % len(c)]][0] * xy[c[i]][1] for i in range(len(c)))))
            for _, c in surf]
    write(base + ".vtk", [c for _, c in surf], "srp_M2.5 mesh (metres)",
          [("quality_minSICN", "double", [f"{q:.5f}" for q in quality]),
           ("cell_size", "double", [f"{s:.6e}" for s in size]),
           ("is_quad", "int", [int(len(c) == 4) for _, c in surf])])

    edges, ids, names = [], [], []
    for i, (name, curve_tags) in enumerate(MARKERS.items()):
        names.append(f"{i}={name}")
        for c in curve_tags:
            for _, e in cells_of(1, c):
                edges.append(e); ids.append(i)
    write(base + "_boundaries.vtk", edges, "srp_M2.5 boundary markers: " + " ".join(names),
          [("marker_id", "int", ids)])


write_vtk_files(OUT.replace(".su2", ""))
gmsh.finalize()
print("wrote", OUT)
