#!/bin/bash
# Mach and density of coarse s6, medium m3, fine f3, fineJet j3, fineJet2 j3 on one 0.025 mm grid over the refinement box,
# for the difference-contour comparison pages (user, 2026-10-30). Run on the desktop node.
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs; P=SA_p0j8.015_p0inf3.247
cd $R/grid_compare || exit 1
unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
H_MM=0.025 pvbatch --mesa resample_box.py box_5meshes_0.025mm.npz coarse=$R/${P}_meshCoarse/2_jet_staged/s6_cap3/flow.vtu \
  medium=$R/${P}_meshMedium/2_jet_staged/m3_cap3/flow.vtu fine=$R/${P}_meshFine/2_jet_staged/f3_cap3/flow.vtu \
  fineJet=$R/${P}_meshFineJet/2_jet_staged/j3_cap3/flow.vtu fineJet2=$R/${P}_meshFineJet2/2_jet_staged/j3_cap3/flow.vtu 2>&1 | grep -vi warning
ls -l box_5meshes_0.025mm.npz
