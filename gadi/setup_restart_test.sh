#!/bin/bash
# Gadi check (2026-11-04): rerun the first N iterations of meshFineJet2 t3 (SST, p0j 7.4945 bar) from a copy of the
# t2 restart, so the history can be compared row by row with M3's t3 and the speed measured. Only the copy of the
# config changes: ITER and the solution file. Usage: bash setup_restart_test.sh [N=500]
N=${1:-500}
R=/scratch/os85/lp3540/su2/cases/srp_M2.5/runs/SST_p0j7.4945_p0inf3.247_meshFineJet2
T=$R/test_gadi
mkdir -p $T && cd $T || exit 1
cp $R/2_jet_staged/t2_cap1/restart.dat start.dat      # a copy, not a symlink (MPI-IO can't read restarts via symlinks)
sed -e "s/^ITER=.*/ITER= $N/" -e "s#^SOLUTION_FILENAME=.*#SOLUTION_FILENAME= start#" \
    $R/2_jet_staged/t3_cap3/t3_cap3.cfg > test.cfg
diff $R/2_jet_staged/t3_cap3/t3_cap3.cfg test.cfg
echo "test folder ready: $T"
