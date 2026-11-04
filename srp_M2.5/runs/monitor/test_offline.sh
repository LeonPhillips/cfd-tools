#!/bin/bash
# Offline test of the monitor (2026-11-04), no SU2: the SST p0j 7.4945 coarse t2 and t3 final solutions stand in for two
# snapshots (symlinks flow_1000.vtu / flow_2000.vtu; the watcher deletes the links, not the solutions).
# Run on a compute node (desktop job): [TEST_NAME=<folder>] bash test_offline.sh. Output: <coarse run>/test_monitor/
# (or TEST_NAME), remove when done.
MON=$(dirname $(readlink -f $0))
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SST_p0j7.4945_p0inf3.247_meshCoarse
T=$R/${TEST_NAME:-test_monitor}; D=$T/2_jet_staged/m_test
[ -e $T ] && { echo "$T exists: remove it first"; exit 1; }
mkdir -p $D
ln -s $R/2_jet_staged/t2_cap1/flow.vtu $D/flow_1000.vtu
ln -s $R/2_jet_staged/t3_cap3/flow.vtu $D/flow_2000.vtu
cp $R/2_jet_staged/t3_cap3/history.csv $D/
( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
  pvbatch --mesa $MON/probes_md.py $R/2_jet_staged/t2_cap1/flow.vtu 2> $D/probes.log | grep -v -i warning | tail -1 > $D/probes.txt )
cat $D/probes.log
{ grep -v "^CUSTOM_OUTPUTS=" $R/2_jet_staged/t3_cap3/t3_cap3.cfg; echo "CUSTOM_OUTPUTS= '$(cat $D/probes.txt)'"; } > $D/m_test.cfg
grep "^CUSTOM_OUTPUTS=" $D/m_test.cfg
ONCE=1 bash $MON/watch_conv.sh $D
echo "== watch.log"; cat $D/watch.log
for f in box_change line_change features; do echo "== $f.csv"; cat $D/$f.csv; done
ls -l $D $D/snaps $T/*_latest.png; du -sh $D
