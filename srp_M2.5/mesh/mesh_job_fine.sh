#!/bin/bash
# meshFine (user, 2026-10-28): the coarse v4 mesh plus a 0.0177 mm box (option B: medium refined by sqrt 2) at x 0.5..20 mm, r 0..9.5 mm (jet, Mach
# disk, bow shock). Everything else as srp_M2.5_coarse_v4 (17 wall layers, first cell 1.8e-7 m). Run on the desktop
# node, not the login node:
#   tmux new -d -s srp_mesh "bash /scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/mesh/mesh_job_fine.sh"
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
cd /scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/mesh
M=srp_M2.5_meshFine
exec > mesh_job_fine.log 2>&1
hostname; date
# 0.025/sqrt(2) = 0.017678 mm failed in Gmsh ("2 intersections in the 1D mesh" at the wall layers, "Edge not recovered"),
# so the sizes nearest to it are tried in turn and the first that meshes is kept (all within 1 % of 0.017678).
for S in 0.0177 0.0176 0.0178 0.0175; do
  BL_LAYERS=17 MED_BOX="0.5,20,0,9.5,$S" python mesh_coarse.py ../geometry/su2_shortDiff_wFillet_v4.STEP $M.su2 > ${M}_mesh.log 2>&1
  rc=$?; echo "box size $S mm: mesh exit=$rc"
  [ $rc -eq 0 ] && break
done
grep -e "extra refinement" -e "2D elements" -e "nodes:" -e "^BL:" -e minSICN -e Error ${M}_mesh.log
python cell_census.py $M.msh
python find_bad_cells.py $M.msh | head -8
sed -e "s/srp_M2.5_coarse/$M/g" meshcheck.cfg > meshcheck_fine.cfg
SU2_CFD meshcheck_fine.cfg > meshcheck_fine.log 2>&1; echo "meshcheck exit=$?"
grep -A6 "Computing mesh quality" meshcheck_fine.log | tail -4; grep "correctly oriented" meshcheck_fine.log
python plot_mesh.py $M.msh ${M}_views.png | tail -1
ls -lh $M*
date
echo MESH_JOB_DONE
