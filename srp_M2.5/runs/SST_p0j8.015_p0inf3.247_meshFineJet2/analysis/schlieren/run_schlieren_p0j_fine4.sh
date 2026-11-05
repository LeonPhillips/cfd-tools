#!/bin/bash
# SST meshFine t3 at p0j 8.015 / 7.4945 / 7.0 / 6.6 bar against the NPR 501 schlieren (user, 2026-11-06).
# 1. Synthetic schlieren for SST_p0j6.6 meshFine t3 (run_schlieren_case.sh; skipped if it exists).
# 2. schlieren_compare_n.py on the four meshFine solutions: one figure in a single column plus one per case.
# Run on a compute/desktop node: srun --jobid=<job> --overlap -n1 -c4 bash <clone>/.../run_schlieren_p0j_fine4.sh
T=$(cd "$(dirname "$0")" && pwd)
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
N=$R/SST_p0j6.6_p0inf3.247_meshFine
[ -f $N/analysis/schlieren/schlieren_SST_p0j6.6_Fine_t3_vs_exp_NPR501.npz ] ||
  bash $T/run_schlieren_case.sh $N 2_jet_staged/t3_cap3 schlieren_SST_p0j6.6_Fine_t3_vs_exp_NPR501 \
    "SST meshFine t3, p0j 6.6 bar" "p0j 6.6 bar, p0inf 3.247 bar" || exit 1
O=$N/analysis/schlieren/four; mkdir -p $O && cd $O || exit 1
module load anaconda/2024.02-1 >/dev/null 2>&1
python $T/schlieren_compare_n.py schlieren_SST_Fine_p0j8.015-7.4945-7.0-6.6_vs_exp_NPR501 \
  "p0j8.015=$R/SST_p0j8.015_p0inf3.247_meshFine/analysis/schlieren/schlieren_SST_Fine_vs_exp_NPR501.npz=SST meshFine t3, p0j 8.015 bar (old inlet)" \
  "p0j7.4945=$R/SST_p0j7.4945_p0inf3.247_meshFine/analysis/schlieren/schlieren_SST_p0j7.4945_Fine_t3_vs_exp_NPR501.npz=SST meshFine t3, p0j 7.4945 bar" \
  "p0j7.0=$R/SST_p0j7.0_p0inf3.247_meshFine/analysis/schlieren/schlieren_SST_p0j7.0_Fine_t3_vs_exp_NPR501.npz=SST meshFine t3, p0j 7.0 bar" \
  "p0j6.6=$N/analysis/schlieren/schlieren_SST_p0j6.6_Fine_t3_vs_exp_NPR501.npz=SST meshFine t3, p0j 6.6 bar"
echo "cfd-tools $(git -C $T rev-parse --short HEAD)" > source.txt
ls -l $O
