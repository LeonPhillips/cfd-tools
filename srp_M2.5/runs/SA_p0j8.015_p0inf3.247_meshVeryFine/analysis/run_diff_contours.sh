#!/bin/bash
# Fine f3 vs very fine v3 difference contours over the refinement box (user, 2026-10-29). Run on the desktop node.
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
F=$R/SA_p0j8.015_p0inf3.247_meshFine; V=$R/SA_p0j8.015_p0inf3.247_meshVeryFine; A=$V/analysis
cd $A || exit 1
unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
pvbatch --mesa resample_box.py box_fine_veryFine.npz fine=$F/2_jet_staged/f3_cap3/flow.vtu veryFine=$V/2_jet_staged/v3_cap3/flow.vtu 2>&1 | grep -vi warning
module unload paraview/5.13 >/dev/null 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null 2>&1
python diff_contours.py box_fine_veryFine.npz $A
ls -l $A/box_fine_veryFine.npz $A/diff_*
