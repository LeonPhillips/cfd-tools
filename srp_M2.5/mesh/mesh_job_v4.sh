#!/bin/bash
# Mesh v4: new geometry su2_shortDiff_wFillet_v4.STEP (outlet at x = -308.1 mm, 0.5 mm base fillet), wall
# layers and sizes as v3, 1.5 mm near-wake box to the outlet (user, 2026-10-27). Run on the desktop node in tmux, not on the login node:
#   tmux new -d -s srp_mesh "bash /scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/mesh/mesh_job_v4.sh"
# v2 failed with 20 wall layers and meshed with 19, so start at 19. On v4, 18 and 19 layers fail ("2 intersections"
# in the wall layers when all walls are one group); 17 meshes, same first cell, layer 0.019 mm instead of 0.028 mm.
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
cd /scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/mesh
M=srp_M2.5_coarse_v4
exec > mesh_job_v4.log 2>&1
hostname; date
for L in 19 18 17; do
  BL_LAYERS=$L python mesh_coarse.py ../geometry/su2_shortDiff_wFillet_v4.STEP $M.su2 > ${M}_mesh.log 2>&1
  rc=$?; echo "BL_LAYERS=$L mesh exit=$rc"
  [ $rc -eq 0 ] && break
done
grep -e "body wall" -e "nozzle throat" -e "used" -e "2D elements" -e "nodes:" -e "^BL:" -e minSICN -e Error ${M}_mesh.log
python cell_census.py $M.msh
python find_bad_cells.py $M.msh | head -8
sed -e "s/srp_M2.5_coarse/$M/g" meshcheck.cfg > meshcheck_v4.cfg
SU2_CFD meshcheck_v4.cfg > meshcheck_v4.log 2>&1; echo "meshcheck exit=$?"
grep -A6 "Computing mesh quality" meshcheck_v4.log | tail -4; grep "correctly oriented" meshcheck_v4.log
python plot_mesh.py $M.msh ${M}_views.png | tail -1
ls -lh $M*
date
echo MESH_JOB_DONE
