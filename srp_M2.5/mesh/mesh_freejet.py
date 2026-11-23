"""Free-jet mesh for the srp_M2.5 nozzle and model, without the counterflow (user, 2026-11-23).

The model and nozzle walls are those of the v4 STEP (curves 1-3 and 10-19, unchanged). The wind-tunnel boundaries
(freestream inlet, slip walls, duct outlet) are replaced by an open rectangular domain:
    x from -308.1241 mm (the v4 rear boundary) to X_OUT (default 240 mm, ~50 exit diameters), r from 0 to R_OUT
    (default 120 mm, ~25 exit diameters).
Markers: axis, inlet_jet (plenum face, as before), wall_noslip (all model walls), ambient_side (r = R_OUT),
ambient_rear (x = -308.1 mm), outlet (x = X_OUT). Found by position after the boolean cut.

Sizes (env, mm): H_NEAR in the near-field box x 0.5..20, r 0..4.5; H_PLUME in the plume boxes x 20..100, whose
radius steps 6.5 / 8.5 / 10.5 / 12 mm every 20 mm (the jet spreads). Levels (user, 2026-11-23, constant sqrt 2 ratio):
    coarse 0.035 / 0.07, medium 0.025 / 0.05, fine 0.0177 / 0.035 (fine = the proposed mesh).
Everything else as meshFine: wall layers (first cell 1.8e-7 m, ratio 1.2, 17 layers), 0.02 mm at the top of the
layers, 0.01 mm on fillets, nozzle + plenum 0.05 mm, 0.6 mm around the body, 1.5 mm behind it, 10 mm far away.
The far jet, x 100..X_OUT, r 0..30 mm, is 0.5 mm on every level (outside the region of interest).
Usage: [H_NEAR=..] [H_PLUME=..] python mesh_freejet.py <v4 step file> <output .su2>
"""
import math
import os
import sys

import gmsh

STEP, OUT = sys.argv[1], sys.argv[2]
mm = 1e-3
X_OUT = float(os.environ.get("X_OUT", "240")) * mm
R_OUT = float(os.environ.get("R_OUT", "120")) * mm
X_REAR = -308.1241 * mm
H_NEAR = float(os.environ.get("H_NEAR", "0.0177")) * mm
H_PLUME = float(os.environ.get("H_PLUME", "0.035")) * mm
H1, BL_RATIO = 1.8e-7, 1.2
BL_LAYERS = int(os.environ.get("BL_LAYERS", "17"))
BL_THICK = H1 * (BL_RATIO ** BL_LAYERS - 1) / (BL_RATIO - 1)
S_FILLET, S_WALL, S_NOZZLE = 1.0e-5, 2.0e-5, 5.0e-5
S_BODY, S_WAKE, S_FARJET, S_FAR = 6.0e-4, 1.5e-3, 5.0e-4, 1.0e-2
GROWTH, REFINE_GROWTH, WAKE_GROWTH = 0.2, 0.1, 0.05
print(f"free jet: x {X_REAR/mm:g}..{X_OUT/mm:g} mm, r 0..{R_OUT/mm:g} mm; near {H_NEAR/mm:g} mm, plume {H_PLUME/mm:g} mm")

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setString("Geometry.OCCTargetUnit", "M")
gmsh.model.add(os.path.basename(OUT).replace(".su2", ""))
gmsh.model.occ.importShapes(STEP)
gmsh.model.occ.synchronize()
assert sorted(t for _, t in gmsh.model.getEntities(1)) == list(range(1, 20)), "not the v4 STEP"


def ends(c):
    return [t for _, t in gmsh.model.getBoundary([(1, c)], oriented=True)]


# body (solid model with the plenum cavity open at the jet inlet face) = model walls + the axis under the body
p_inlet_axis, p_base_axis = ends(3)[1], ends(10)[0]          # (-8.1241, 0) and (-58.1241, 0)
old = [t for _, t in gmsh.model.getEntities(2)]
gmsh.model.occ.remove([(2, t) for t in old])
under = gmsh.model.occ.addLine(p_inlet_axis, p_base_axis)
loop = gmsh.model.occ.addCurveLoop([3, under, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 1, 2])
body = gmsh.model.occ.addPlaneSurface([loop])
box_dom = gmsh.model.occ.addRectangle(X_REAR, 0, 0, X_OUT - X_REAR, R_OUT)
fluid, _ = gmsh.model.occ.cut([(2, box_dom)], [(2, body)])
gmsh.model.occ.synchronize()
# drop leftovers (the old curves not on the new fluid boundary)
keep = {t for _, t in gmsh.model.getBoundary(fluid, oriented=False)}
gmsh.model.occ.remove([(1, t) for _, t in gmsh.model.getEntities(1) if t not in keep], recursive=True)
gmsh.model.occ.synchronize()
assert len(gmsh.model.getEntities(2)) == 1, gmsh.model.getEntities(2)
surf = gmsh.model.getEntities(2)[0][1]

tol = 1e-6   # bounding boxes are padded by 1e-7 m
MARKERS = {k: [] for k in ("axis", "inlet_jet", "wall_noslip", "ambient_side", "ambient_rear", "outlet")}
fillets = []
for _, c in gmsh.model.getBoundary([(2, surf)], oriented=False):
    x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(1, c)
    if abs(y0) < tol and abs(y1) < tol:
        MARKERS["axis"].append(c)
    elif abs(y0 - R_OUT) < tol and abs(y1 - R_OUT) < tol:
        MARKERS["ambient_side"].append(c)
    elif abs(x0 - X_REAR) < tol and abs(x1 - X_REAR) < tol:
        MARKERS["ambient_rear"].append(c)
    elif abs(x0 - X_OUT) < tol and abs(x1 - X_OUT) < tol:
        MARKERS["outlet"].append(c)
    elif abs(x0 + 8.1241 * mm) < 1e-6 and abs(x1 + 8.1241 * mm) < 1e-6 and y1 < 2.9 * mm:
        MARKERS["inlet_jet"].append(c)
    else:
        MARKERS["wall_noslip"].append(c)
        if gmsh.model.getType(1, c) == "Circle":
            fillets.append(c)
for k, v in MARKERS.items():
    xs = [gmsh.model.getBoundingBox(1, c) for c in v]
    print(f"{k:13s} curves {sorted(v)}  x {min(b[0] for b in xs)/mm:.4f}..{max(b[3] for b in xs)/mm:.4f}  "
          f"r {min(b[1] for b in xs)/mm:.4f}..{max(b[4] for b in xs)/mm:.4f} mm")
    assert v, f"no curves for {k}"
assert len(MARKERS["inlet_jet"]) == 1 and len(MARKERS["wall_noslip"]) >= 12, "unexpected marker split"
for name, tags in MARKERS.items():
    gmsh.model.addPhysicalGroup(1, tags, name=name)
gmsh.model.addPhysicalGroup(2, [surf], name="fluid")
if os.environ.get("GEOM_ONLY") == "1":   # quick check of the geometry and markers, no mesh
    gmsh.finalize(); sys.exit(0)

f = gmsh.model.mesh.field
walls = MARKERS["wall_noslip"]


def point_at(x, y):  # nearest geometry point, which must be within 1 micron
    d, p = min((math.hypot(gmsh.model.getValue(0, p, [])[0] - x, gmsh.model.getValue(0, p, [])[1] - y), p)
               for _, p in gmsh.model.getEntities(0))
    assert d < 1e-6, f"no point near {x}, {y} (nearest {d} m)"
    return p


# wall layers end cleanly where they meet the axis (base) and the jet inlet (plenum wall), as in mesh_coarse.py
bl_end = [point_at(-58.1241 * mm, 0.0), point_at(-8.1241 * mm, 2.8787 * mm)]
bl = f.add("BoundaryLayer")
f.setNumbers(bl, "CurvesList", walls)
f.setNumber(bl, "Size", H1)
f.setNumber(bl, "Ratio", BL_RATIO)
f.setNumber(bl, "Thickness", BL_THICK)
f.setNumber(bl, "Quads", 1)
f.setNumbers(bl, "PointsList", bl_end)
f.setAsBoundaryLayer(bl)


def threshold(curve_list, s_min, s_max, d_min, d_max, local=False):
    d = f.add("Distance")
    f.setNumbers(d, "CurvesList", curve_list)
    f.setNumber(d, "Sampling", 400)
    t = f.add("Threshold")
    f.setNumber(t, "InField", d)
    for k, v in dict(SizeMin=s_min, SizeMax=s_max, DistMin=d_min, DistMax=d_max).items():
        f.setNumber(t, k, v)
    f.setNumber(t, "StopAtDistMax", 1 if local else 0)
    return t


def box(x0, x1, y0, y1, size, v_out=S_FAR, growth=GROWTH):
    b = f.add("Box")
    for k, v in dict(XMin=x0, XMax=x1, YMin=y0, YMax=y1, ZMin=-1, ZMax=1,
                     VIn=size, VOut=v_out, Thickness=(v_out - size) / growth).items():
        f.setNumber(b, k, v)
    return b


fields = [
    threshold(fillets, S_FILLET, S_WALL, 2e-5, 2e-4, local=True),
    threshold(walls, S_WALL, S_FAR, BL_THICK, BL_THICK + (S_FAR - S_WALL) / GROWTH),
    box(-8.2 * mm, 2.0 * mm, 0.0, 3.0 * mm, S_NOZZLE),                       # nozzle + plenum
    box(0.5 * mm, 20 * mm, 0.0, 4.5 * mm, H_NEAR, growth=REFINE_GROWTH),      # near field: lip, shear layer, cells
    box(100 * mm, X_OUT, 0.0, 30 * mm, S_FARJET),                             # far jet
    box(-60 * mm, 100 * mm, 0.0, 30 * mm, S_BODY),                            # around the body and the jet
    box(X_REAR - 2 * mm, -10 * mm, 0.0, 60 * mm, S_WAKE, growth=WAKE_GROWTH),  # behind the body
]
for i, r in enumerate((6.5, 8.5, 10.5, 12.0)):                               # plume, widening every 20 mm
    fields.append(box((20 + 20 * i) * mm, (40 + 20 * i) * mm, 0.0, r * mm, H_PLUME, growth=REFINE_GROWTH))
fmin = f.add("Min")
f.setNumbers(fmin, "FieldsList", fields)
f.setAsBackgroundMesh(fmin)

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.Algorithm", int(os.environ.get("MESH_ALGO", "8")))
gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 0)
gmsh.option.setNumber("Mesh.Smoothing", 5)
gmsh.model.mesh.generate(2)

types, tags, _ = gmsh.model.mesh.getElements(2)
counts = {gmsh.model.mesh.getElementProperties(t)[0]: len(e) for t, e in zip(types, tags)}
print("2D elements:", counts, "total", sum(counts.values()))
print("nodes:", len(gmsh.model.mesh.getNodes()[0]))
q = gmsh.model.mesh.getElementQualities([e for tl in tags for e in tl], "minSICN")
print(f"minSICN: min={min(q):.3f}  <0.1: {sum(x < 0.1 for x in q)}  <0: {sum(x <= 0 for x in q)}")
print(f"BL: first cell {H1:.3e} m, {BL_LAYERS} layers, ratio {BL_RATIO}, thickness {BL_THICK:.3e} m")
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(OUT)
gmsh.write(OUT.replace(".su2", ".msh"))

# marker map figure data: boundary edges as legacy VTK (marker_id), as mesh_coarse.py
ntags, coords, _ = gmsh.model.mesh.getNodes()
index = {t: i for i, t in enumerate(ntags)}
xyz = coords.reshape(-1, 3)
edges, ids, names = [], [], []
for i, (name, ctags) in enumerate(MARKERS.items()):
    names.append(f"{i}={name}")
    for c in ctags:
        for etype, _, en in zip(*gmsh.model.mesh.getElements(1, c)):
            n = gmsh.model.mesh.getElementProperties(etype)[3]
            for e in en.reshape(-1, n):
                edges.append([index[k] for k in e]); ids.append(i)
with open(OUT.replace(".su2", "_boundaries.vtk"), "w") as fh:
    fh.write(f"# vtk DataFile Version 3.0\nsrp_M2.5 free-jet boundary markers: {' '.join(names)}\nASCII\n"
             f"DATASET UNSTRUCTURED_GRID\nPOINTS {len(xyz)} double\n")
    fh.write("\n".join(f"{x:.9e} {y:.9e} 0" for x, y, _ in xyz) + "\n")
    fh.write(f"CELLS {len(edges)} {3 * len(edges)}\n" + "\n".join(f"2 {a} {b}" for a, b in edges) + "\n")
    fh.write(f"CELL_TYPES {len(edges)}\n" + "\n".join("3" for _ in edges) + "\n")
    fh.write(f"CELL_DATA {len(edges)}\nSCALARS marker_id int 1\nLOOKUP_TABLE default\n" + "\n".join(map(str, ids)) + "\n")
gmsh.finalize()
print("wrote", OUT)
