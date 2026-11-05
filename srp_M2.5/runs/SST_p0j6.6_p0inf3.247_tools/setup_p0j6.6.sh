#!/bin/bash
# Build SST_p0j6.6_p0inf3.247_* (user, 2026-11-05): the SST p0j 7.0 bar study repeated at jet total pressure
# p0j = 6.6 bar (T0j 288.15 K, TOTAL_CONDITIONS), on meshes coarse, medium, fine and veryFine.
# Configs are the SST_p0j7.0_p0inf3.247 ones (all four meshes exist there) with only the inlet line, headers and
# sources changed.
# Starts: coarse restarts from the 7.0 bar coarse final solution (t3_cap3, same mesh). Every other mesh starts from
# the previous, coarser mesh's new 6.6 bar t3, interpolated (0_interp): medium <- coarse, fine <- medium,
# veryFine <- fine. Then t1 cap 0.5 (5000) / t2 cap 1 (25000) / t3 cap 3 (20000). Nothing stops on convergence.
# The probe line (CUSTOM_OUTPUTS) is set by the chain before each stage (monitor/probes_md.py).
# Run on M3 from the cfd-tools clone: bash <clone>/srp_M2.5/runs/SST_p0j6.6_p0inf3.247_tools/setup_p0j6.6.sh
set -e
SRC=$(dirname $(readlink -f $0))
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
OLD=SST_p0j7.0_p0inf3.247; NEW=SST_p0j6.6_p0inf3.247
MESHES="Coarse Medium Fine VeryFine"
declare -A FROM=([Medium]=Coarse [Fine]=Medium [VeryFine]=Fine)
for m in $MESHES; do [ -e $R/${NEW}_mesh$m ] && { echo "$R/${NEW}_mesh$m exists: not touching it"; exit 1; }; done
[ -e $R/${NEW}_tools ] && { echo "$R/${NEW}_tools exists"; exit 1; }
[ -s $R/${OLD}_meshCoarse/2_jet_staged/t3_cap3/restart.dat ] || { echo "missing ${OLD}_meshCoarse t3_cap3/restart.dat"; exit 1; }

# tools: this folder plus the monitor and the feature finder (same relative layout, so the monitor's imports work)
mkdir -p $R/${NEW}_tools/monitor $R/${NEW}_tools/grid_compare
cp $SRC/*.sh $SRC/README.md $R/${NEW}_tools/
cp $SRC/../monitor/*.py $SRC/../monitor/*.sh $SRC/../monitor/README.md $R/${NEW}_tools/monitor/
cp $SRC/../grid_compare/feature_positions.py $R/${NEW}_tools/grid_compare/
echo "copied from cfd-tools $(git -C $SRC rev-parse HEAD) on $(date)" > $R/${NEW}_tools/COMMIT

INLET='MARKER_INLET= ( inlet_jet, 288.15, 660000.0, 1.0, 0.0, 0.0 )   % p0j 6.6 bar (total), T0j 288.15 K (user, 2026-11-05)'
NOTE='% Phase 2: inlet_jet is a total-conditions inlet, p0j = 6.6 bar, T0j = 288.15 K, flow in +x (user 2026-11-05).'
edit() {  # edit <header line> <source cfg>: the common changes, to stdout
  sed -e "s|$OLD|$NEW|g" \
      -e "s|^% Nothing stops on convergence. Made by.*|% Nothing stops on convergence. Made by ../../../${NEW}_tools/setup_p0j6.6.sh from the $OLD configs.|" \
      -e "s|^% Phase 2: inlet_jet is.*|$NOTE|" \
      -e "s|^MARKER_INLET=.*|$INLET|" \
      -e "1s|.*|$1|" $2
}
cd $R
for m in $MESHES; do
  t=${OLD}_mesh$m
  mesh=$(sed -n 's/^MESH_FILENAME= *//p' $t/2_jet_staged/t3_cap3/t3_cap3.cfg)
  [ -s $mesh ] || { echo "missing mesh $mesh"; exit 1; }
  if [ $m = Coarse ]; then
    prev="from $OLD meshCoarse t3_cap3 (7.0 bar inlet), same mesh"
  else
    g=${NEW}_mesh$m/0_interp/interp.cfg; mkdir -p $(dirname $g)
    edit "% srp_M2.5, case $NEW, mesh$m, k-omega SST (user, 2026-11-05). Interpolation of mesh${FROM[$m]} t3_cap3 (1 iteration at CFL 1e-8)" \
      $t/0_interp/interp.cfg > $g
    sed -i "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $R/${NEW}_mesh${FROM[$m]}/2_jet_staged/t3_cap3/restart|" $g
    prev="from mesh${FROM[$m]} t3_cap3, interpolated (0_interp)"
  fi
  for s in t1_cap0.5 t2_cap1 t3_cap3; do
    g=${NEW}_mesh$m/2_jet_staged/$s/$s.cfg; mkdir -p $(dirname $g)
    edit "% srp_M2.5, case $NEW, mesh$m, k-omega SST (user, 2026-11-05). Stage $s: $prev" $t/2_jet_staged/$s/$s.cfg > $g
    if [ $s = t1_cap0.5 ]; then
      [ $m = Coarse ] && src=$R/${OLD}_meshCoarse/2_jet_staged/t3_cap3/restart || src=../../0_interp/restart
      sed -i "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $src|" $g
    fi
    prev="from $s"
  done
done
CF=$(find ${NEW}_mesh* -name "*.cfg" | sort)
echo "built:"; echo "$CF"
echo "inlet lines:"; grep -h "^MARKER_INLET=\|^INLET_TYPE" $CF | sort | uniq -c
echo "mesh files:"; for m in $MESHES; do echo "  $m $(grep -h '^MESH_FILENAME' $(find ${NEW}_mesh$m -name '*.cfg') | sort -u)"; done
echo "sources:"; grep -H "^SOLUTION_FILENAME" $CF
echo "ITER / output:"; grep -H "^ITER=\|^OUTPUT_WRT_FREQ=\|^WRT_VOLUME_OVERWRITE=" $CF | sed 's/ *%.*//'
echo "old-case or 7.0 bar references left (header line 3 and the coarse source are expected):"
grep -n "$OLD\|700000\|7\.0 bar" $CF | grep -v "^[^:]*:3:\|SOLUTION_FILENAME= $R/${OLD}_meshCoarse\|^[^:]*:1:.*7.0 bar inlet" || echo none
