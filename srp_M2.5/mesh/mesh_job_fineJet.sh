#!/bin/bash
# meshFineJet (user, 2026-10-29): meshFine (0.0177 mm box at x 0.5..20 mm, r 0..9.5 mm) plus a local 0.0125 mm box at
# x 1.5..8.5 mm, r 0..3.6 mm around the lip, shear layer, Mach disk and terminal shock, where fine and very fine still
# differ. Nozzle stays at 0.05 mm. 17 wall layers, first cell 1.8e-7 m. Run on the desktop node.
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
cd /scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/mesh
M=srp_M2.5_meshFineJet
exec > mesh_job_fineJet.log 2>&1
hostname; date
# Fallback sizes within 1 % in case Gmsh hits the "2 intersections" wall-layer failure.
for S in 0.0125 0.0124 0.0126; do
  BL_LAYERS=17 MED_BOX="0.5,20,0,9.5,0.0177;1.5,8.5,0,3.6,$S" python mesh_coarse.py ../geometry/su2_shortDiff_wFillet_v4.STEP $M.su2 > ${M}_mesh.log 2>&1
  rc=$?; echo "jet box size $S mm: mesh exit=$rc"
  [ $rc -eq 0 ] && break
done
grep -e "extra refinement" -e "2D elements" -e "nodes:" -e "^BL:" -e minSICN -e Error ${M}_mesh.log
python cell_census.py $M.msh
python find_bad_cells.py $M.msh | head -8
sed -e "s/srp_M2.5_coarse/$M/g" meshcheck.cfg > meshcheck_fineJet.cfg
SU2_CFD meshcheck_fineJet.cfg > meshcheck_fineJet.log 2>&1; echo "meshcheck exit=$?"
grep -A6 "Computing mesh quality" meshcheck_fineJet.log | tail -4; grep "correctly oriented" meshcheck_fineJet.log
python plot_mesh.py $M.msh ${M}_views.png | tail -1
python compare_fine_fineJet.py fine_vs_fineJet.png
ls -lh $M* fine_vs_fineJet.png
date
echo MESH_JOB_DONE
