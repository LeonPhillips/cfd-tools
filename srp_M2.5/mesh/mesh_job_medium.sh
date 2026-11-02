#!/bin/bash
# meshMedium (user, 2026-10-28): the coarse v4 mesh plus a 0.025 mm box at x 0.5..20 mm, r 0..9.5 mm (jet, Mach
# disk, bow shock). Everything else as srp_M2.5_coarse_v4 (17 wall layers, first cell 1.8e-7 m). Run on the desktop
# node, not the login node:
#   tmux new -d -s srp_mesh "bash /scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/mesh/mesh_job_medium.sh"
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
cd /scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/mesh
M=srp_M2.5_meshMedium
exec > mesh_job_medium.log 2>&1
hostname; date
BL_LAYERS=17 MED_BOX="0.5,20,0,9.5,0.025" python mesh_coarse.py ../geometry/su2_shortDiff_wFillet_v4.STEP $M.su2 > ${M}_mesh.log 2>&1
echo "mesh exit=$?"
grep -e "extra refinement" -e "2D elements" -e "nodes:" -e "^BL:" -e minSICN -e Error ${M}_mesh.log
python cell_census.py $M.msh
python find_bad_cells.py $M.msh | head -8
sed -e "s/srp_M2.5_coarse/$M/g" meshcheck.cfg > meshcheck_medium.cfg
SU2_CFD meshcheck_medium.cfg > meshcheck_medium.log 2>&1; echo "meshcheck exit=$?"
grep -A6 "Computing mesh quality" meshcheck_medium.log | tail -4; grep "correctly oriented" meshcheck_medium.log
python plot_mesh.py $M.msh ${M}_views.png | tail -1
ls -lh $M*
date
echo MESH_JOB_DONE
