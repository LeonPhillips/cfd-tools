#!/bin/bash
# Synthetic d(rho)/dx schlieren from SST meshFineJet2 t3 vs the mean experimental image (user, 2026-11-02).
# Run on the desktop node. Experimental extract: x, y (in D), im_mean, background from the Market file
# SCH1_ZOOM1_TOP_CTa_04891_NPRa_463_CTd_04925_NPRd_467_CT_04530_NPR_501_220_800_1_Fields.mat
A=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SST_p0j8.015_p0inf3.247_meshFineJet2/analysis/schlieren
cd $A || exit 1
( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
  [ -f cfd_drhodx.npz ] || pvbatch --mesa sample_schlieren.py ../../2_jet_staged/t3_cap3/flow.vtu cfd_drhodx.npz 2>&1 | grep -vi warning )
module load anaconda/2024.02-1 >/dev/null 2>&1
python schlieren_compare.py cfd_drhodx.npz SCH1_ZOOM1_NPR501_220_800_1_mean.mat schlieren_SST_fineJet2_vs_exp_NPR501
ls -l $A
