#!/bin/bash
# meshVeryFine (user, 2026-10-28): the coarse v4 mesh plus a 0.0125 mm box (option A: fine refined by sqrt 2, so medium /
# fine / very fine are 0.025 / 0.0177 / 0.0125 mm, a constant ratio) at x 0.5..20 mm, r 0..9.5 mm. The nozzle stays at
# 0.05 mm (meshFineNozzle showed that is adequate). 17 wall layers, first cell 1.8e-7 m. Run on the desktop node.
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
cd /scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/mesh
M=srp_M2.5_meshVeryFine
exec > mesh_job_veryFine.log 2>&1
hostname; date
# Fallback sizes within 1 % in case Gmsh hits the "2 intersections" wall-layer failure (0.0125 meshed on 2026-10-28).
for S in 0.0125 0.0124 0.0126; do
  BL_LAYERS=17 MED_BOX="0.5,20,0,9.5,$S" python mesh_coarse.py ../geometry/su2_shortDiff_wFillet_v4.STEP $M.su2 > ${M}_mesh.log 2>&1
  rc=$?; echo "box size $S mm: mesh exit=$rc"
  [ $rc -eq 0 ] && break
done
grep -e "extra refinement" -e "2D elements" -e "nodes:" -e "^BL:" -e minSICN -e Error ${M}_mesh.log
python cell_census.py $M.msh
python find_bad_cells.py $M.msh | head -8
sed -e "s/srp_M2.5_coarse/$M/g" meshcheck.cfg > meshcheck_veryFine.cfg
SU2_CFD meshcheck_veryFine.cfg > meshcheck_veryFine.log 2>&1; echo "meshcheck exit=$?"
grep -A6 "Computing mesh quality" meshcheck_veryFine.log | tail -4; grep "correctly oriented" meshcheck_veryFine.log
python plot_mesh.py $M.msh ${M}_views.png | tail -1
python compare_fine_veryFine.py fine_vs_veryFine.png
ls -lh $M* fine_vs_veryFine.png
date
echo MESH_JOB_DONE
