#!/bin/bash
# meshFineNozzle (user, 2026-10-28): meshFine (0.0177 mm box at x 0.5..20 mm, r 0..9.5 mm) plus the nozzle + plenum
# interior x -8.2..0.5 mm, r 0..3 mm refined from 0.05 to the same 0.0177 mm (option B). 17 wall layers, first cell 1.8e-7 m.
# Run on the desktop node, not the login node.
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
cd /scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/mesh
M=srp_M2.5_meshFineNozzle
exec > mesh_job_fineNozzle.log 2>&1
hostname; date
# Same fallback as meshFine: sizes within 1 % of 0.0177 mm in case Gmsh hits the "2 intersections" wall-layer failure.
for S in 0.0177 0.0176 0.0178 0.0175; do
  BL_LAYERS=17 MED_BOX="0.5,20,0,9.5,$S;-8.2,0.5,0,3,$S" python mesh_coarse.py ../geometry/su2_shortDiff_wFillet_v4.STEP $M.su2 > ${M}_mesh.log 2>&1
  rc=$?; echo "box size $S mm: mesh exit=$rc"
  [ $rc -eq 0 ] && break
done
grep -e "extra refinement" -e "2D elements" -e "nodes:" -e "^BL:" -e minSICN -e Error ${M}_mesh.log
python cell_census.py $M.msh
python find_bad_cells.py $M.msh | head -8
sed -e "s/srp_M2.5_coarse/$M/g" meshcheck.cfg > meshcheck_fineNozzle.cfg
SU2_CFD meshcheck_fineNozzle.cfg > meshcheck_fineNozzle.log 2>&1; echo "meshcheck exit=$?"
grep -A6 "Computing mesh quality" meshcheck_fineNozzle.log | tail -4; grep "correctly oriented" meshcheck_fineNozzle.log
python plot_mesh.py $M.msh ${M}_views.png | tail -1
python compare_fine_nozzle.py fine_vs_fineNozzle.png
ls -lh $M* fine_vs_fineNozzle.png
date
echo MESH_JOB_DONE
