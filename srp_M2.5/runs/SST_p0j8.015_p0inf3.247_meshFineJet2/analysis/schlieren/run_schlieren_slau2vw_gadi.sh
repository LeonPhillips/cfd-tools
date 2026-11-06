#!/bin/bash
# SST meshFine SLAU2 + Venkatakrishnan-Wang (end of t5, cumulative iteration 99999) at p0j 7.4945 / 7.0 / 6.6 bar against
# the NPR 501 schlieren, on Gadi (user, 2026-11-06). The same steps as run_schlieren_case.sh + run_schlieren_p0j_fine4.sh
# on M3, with Gadi paths and ParaView 6.1.1 (pvbatch has no --mesa flag; it runs in a subshell without the SU2 environment).
# 1. per case: sample_schlieren.py (d(rho)/dx, Gradient filter) -> schlieren_compare.py -> <run>/analysis/schlieren/
# 2. schlieren_compare_n.py on the three cases -> SST_p0j6.6_..._SLAU2-VW/analysis/schlieren/three/
# Run inside a job: qsub schlieren_slau2vw.pbs
T=$(cd "$(dirname "$0")" && pwd)
R=/scratch/os85/lp3540/su2/cases/srp_M2.5/runs
EXP=/scratch/os85/lp3540/su2/cases/srp_M2.5/experiment/SCH1_ZOOM1_NPR501_220_800_1_mean.mat
PV_MODULE=${PV_MODULE:-paraview/6.1.1-mesa}
source /scratch/os85/lp3540/su2/SU2-v8.5.0/su2_env.sh >/dev/null
type module >/dev/null 2>&1 || source /opt/Modules/v4.3.0/init/bash
COMMIT=$(git -C $T rev-parse --short HEAD)
args=()
for P in 7.4945 7.0 6.6; do
  RUN=$R/SST_p0j${P}_p0inf3.247_meshFine_SLAU2-VW; STAGE=2_jet_staged/t5_cap3
  A=$RUN/analysis/schlieren; PRE=schlieren_SST_p0j${P}_Fine_SLAU2-VW_t5_vs_exp_NPR501
  mkdir -p $A && cd $A || exit 1
  echo "== p0j $P $(date +%T)"
  [ -f cfd_drhodx.npz ] || ( unset LD_PRELOAD PYTHONPATH PYTHONHOME VIRTUAL_ENV; module purge >/dev/null 2>&1
      module load $PV_MODULE >/dev/null 2>&1; pvbatch $T/sample_schlieren.py $RUN/$STAGE/flow.vtu cfd_drhodx.npz 2>&1 | grep -vi warning )
  [ -f cfd_drhodx.npz ] || { echo "sampling failed for p0j $P"; exit 1; }
  CFD_LABEL="SST meshFine SLAU2-VW t5" CASE_LABEL="p0j $P bar, p0inf 3.247 bar" python $T/schlieren_compare.py cfd_drhodx.npz $EXP $PRE || exit 1
  echo "cfd-tools $COMMIT, stage $STAGE (iteration 49999 of t5, cumulative 99999)" > source.txt
  args+=("p0j$P=$A/$PRE.npz=SST meshFine SLAU2-VW t5, p0j $P bar")
done
O=$R/SST_p0j6.6_p0inf3.247_meshFine_SLAU2-VW/analysis/schlieren/three; mkdir -p $O && cd $O || exit 1
python $T/schlieren_compare_n.py schlieren_SST_Fine_SLAU2-VW_p0j7.4945-7.0-6.6_vs_exp_NPR501 "${args[@]}" || exit 1
echo "cfd-tools $COMMIT" > source.txt
ls -l $O
