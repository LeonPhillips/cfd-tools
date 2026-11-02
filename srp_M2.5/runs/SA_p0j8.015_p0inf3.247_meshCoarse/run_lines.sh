#!/bin/bash
# Line profiles (Mach, density) for the grid comparison and the fine f2-vs-f3 convergence check. Run on the desktop node.
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
C=$R/SA_p0j8.015_p0inf3.247_meshCoarse; MD=$R/SA_p0j8.015_p0inf3.247_meshMedium; F=$R/SA_p0j8.015_p0inf3.247_meshFine
G=$R/grid_compare; mkdir -p $G; cd $G || exit 1
unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
pvbatch --mesa $C/line_profiles.py grid_lines.npz coarse=$C/2_jet_staged/s6_cap3/flow.vtu \
  medium=$MD/2_jet_staged/m3_cap3/flow.vtu fine=$F/2_jet_staged/f3_cap3/flow.vtu 2>&1 | grep -vi warning
pvbatch --mesa $C/line_profiles.py $F/analysis/conv_lines_f2_f3.npz f2_end=$F/2_jet_staged/f2_cap1/flow.vtu \
  f3_end=$F/2_jet_staged/f3_cap3/flow.vtu 2>&1 | grep -vi warning
module unload paraview/5.13 >/dev/null 2>&1
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null 2>&1
echo "== grid comparison"
python $C/plot_lines.py grid_lines.png grid_lines.npz "Final solutions: coarse s6, medium m3, fine f3"
echo "== convergence check: fine, end of f2 vs end of f3 (20000 iterations apart)"
python $C/plot_lines.py $F/analysis/conv_lines_f2_f3.png $F/analysis/conv_lines_f2_f3.npz \
  "meshFine: end of f2 (it 25000) vs end of f3 (it 20000 later)"
ls -l $G $F/analysis/conv_lines*
