#!/bin/bash
# Short checks before the p0j 7.0 bar chain (user, 2026-11-04). Run after setup_p0j7.0.sh, inside the desktop job:
#   srun --jobid=<job> --overlap -n1 -c24 bash <run-tree tools>/test_p0j7.0.sh      (NP ranks, default 22)
# 1. coarse t1 (the built config) for ITER=${ITER:-2000} from the 7.4945 bar coarse t3, with the probes from
#    probes_md.py and the live monitor: checks the inlet line, the jet mass flow (expected ~0.0207 kg/s once the plenum
#    has settled), the probes and the monitor figures.
# 2. veryFine interpolation (the built interp config, with the 7.4945 bar fine t3 standing in for the 7.0 bar fine t3),
#    then probes_md.py on it.
# 3. plot_grid_lines.py / plot_grid_features.py on monitor snapshots of the 7.4945 bar coarse, medium and fine t3 and
#    the interpolated veryFine field from 2 (snap.py + conv_update.py, no SU2).
# Output: runs/SST_p0j7.0_p0inf3.247_test/ (remove when done). The built run folders are not touched.
TOOLS=$(dirname $(readlink -f $0)); MON=$TOOLS/monitor
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
NEW=SST_p0j7.0_p0inf3.247; OLD=SST_p0j7.4945_p0inf3.247; T=$R/${NEW}_test; NP=${NP:-22}
[ -e $T ] && { echo "$T exists: remove it first"; exit 1; }
mkdir -p $T; cd $T
echo "host $(hostname), cpus $(grep Cpus_allowed_list /proc/self/status | cut -f2), $NP ranks, $(date)"
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null
export I_MPI_HYDRA_BOOTSTRAP=fork
pv() { ( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1; pvbatch --mesa "$@" ); }
mflow() { python - $1 <<'PY'
import csv, sys
r = list(csv.reader(open(sys.argv[1]))); h = [c.strip().strip('"') for c in r[0]]
for k in ("Avg_Massflow(inlet_jet)", "Avg_Massflow(inlet_freestream)", "Avg_Massflow(outlet)"):
    if k in h:
        print(f"  {k}: " + ", ".join(f"it {int(float(x[2]))} {float(x[h.index(k)]):.6f}" for x in r[1::max(1, len(r) // 5)][:6] + [r[-1]]))
PY
}

echo "== 1. coarse t1, ITER ${ITER:-2000}"
D=$T/coarse_t1/2_jet_staged/t1_cap0.5; mkdir -p $D; cd $D
sed "s|^ITER=.*|ITER= ${ITER:-2000}|" $R/${NEW}_meshCoarse/2_jet_staged/t1_cap0.5/t1_cap0.5.cfg > t1_cap0.5.cfg
p=$(pv $MON/probes_md.py $R/${OLD}_meshCoarse/2_jet_staged/t3_cap3/flow.vtu 2> probes.log | grep -v -i warning | tail -1)
cat probes.log
[[ "$p" == up0* ]] || { echo "probes_md.py failed"; exit 1; }
sed -i "s|^CUSTOM_OUTPUTS=.*|CUSTOM_OUTPUTS= '$p'|" t1_cap0.5.cfg
grep -E "^(MARKER_INLET|INLET_TYPE|SOLUTION_FILENAME|MESH_FILENAME|ITER|OUTPUT_WRT_FREQ|WRT_VOLUME_OVERWRITE)=" t1_cap0.5.cfg
grep -c "up0 : Probe" t1_cap0.5.cfg | sed 's/^/  CUSTOM_OUTPUTS lines with the new probes: /'
(sleep 30; HOLD=0 PREV=- bash $MON/watch_conv.sh $D) &
{ echo "host $(hostname)  ranks $NP"; date; } > t1_cap0.5.log
mpirun -n $NP SU2_CFD -t 1 t1_cap0.5.cfg >> t1_cap0.5.log 2>&1
echo "SU2_CFD exit=$?  $(date)" | tee -a t1_cap0.5.log
wait
mflow history.csv
echo "  (7.4945 bar run: jet 0.022165 kg/s; expected at 7.0 bar ~0.0207)"
echo "-- watch.log"; cat watch.log
for f in line_change features; do echo "-- $f.csv"; cat $f.csv; done
ls -l $D/*.png $T/coarse_t1/*_latest.png

echo "== 2. veryFine interpolation (7.4945 bar fine t3 standing in)"
V=$T/veryFine_interp; mkdir -p $V; cd $V
sed "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $R/${OLD}_meshFine/2_jet_staged/t3_cap3/restart|" $R/${NEW}_meshVeryFine/0_interp/interp.cfg > interp.cfg
grep -E "^(MARKER_INLET|SOLUTION_FILENAME|MESH_FILENAME|ITER)=" interp.cfg
{ echo "host $(hostname)  ranks $NP"; date; } > interp.log
mpirun -n $NP SU2_CFD -t 1 interp.cfg >> interp.log 2>&1
echo "SU2_CFD exit=$?  $(date)" | tee -a interp.log
grep -q "does not match the mesh" interp.log && echo "  (SU2 interpolated: point counts differ)"
ls -l restart.dat flow.vtu
pv $MON/probes_md.py $V/flow.vtu 2>&1 >/dev/null | grep -v -i warning

echo "== 3. grid plots"
for m in Coarse Medium Fine VeryFine; do
  G=$T/grid/$m; mkdir -p $G/snaps
  [ $m = VeryFine ] && f=$V/flow.vtu || f=$R/${OLD}_mesh$m/2_jet_staged/t3_cap3/flow.vtu
  ln -s $f $G/flow_19999.vtu
  pv $MON/snap.py $G/snaps $G/flow_19999.vtu 2>&1 | grep -vi warning
  rm $G/flow_19999.vtu
  python $MON/conv_update.py $G
done
cd $T/grid
A="coarse=Coarse medium=Medium fine=Fine veryFine-interp=VeryFine"
python $MON/plot_grid_lines.py grid_lines_test.png "TEST: SST p0j 7.4945 bar t3 (coarse-fine) + fine interpolated onto veryFine" $A
python $MON/plot_grid_features.py grid_features_test.png "TEST: SST p0j 7.4945 bar t3 (coarse-fine) + fine interpolated onto veryFine" $A
cat grid_features_test.csv
ls -l $T/grid/*.png; du -sh $T
