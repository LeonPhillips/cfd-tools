#!/bin/bash
# Build SST_p0j7.4945_p0inf3.247_* (user, 2026-11-02): the SST mesh study (SST_p0j8.015_p0inf3.247_*) repeated with the
# corrected jet inlet. The user gave static p = 7.0588 bar and V = 99.11 m/s at the jet inlet (plenum face); with
# T0 = 288.15 K (option A): T = 288.15 - 99.11^2/(2*1004.703) = 283.26 K, M = 0.2937, p0 = 7.0588*(1 + 0.2 M^2)^3.5
# = 7.4945 bar. The inlet stays TOTAL_CONDITIONS; only the p0j value changes.
# Coarse restarts from its own SST 8.015 bar final solution (t4_cap3, same mesh). Every other mesh starts from the new
# final (t3_cap3) solution of the previous, coarser mesh, interpolated (0_interp; user, 2026-11-03): medium <- coarse,
# fine <- medium, fineJet <- fine, fineJet2 <- fineJet, veryFine <- fineJet2 (user, 2026-11-03; the SST study used fine).
# Then t1 cap 0.5 (5000) / t2 cap 1 (25000) / t3 cap 3 (20000) on every mesh. The coarse configs are the meshMedium ones
# with the coarse mesh file. Nothing stops on convergence.
# Run on M3 from the cfd-tools clone: bash <clone>/srp_M2.5/runs/SST_p0j7.4945_p0inf3.247_tools/setup_p0j7.4945.sh
# UPDATE=1 (2026-11-03, while the chain was running): leave coarse and any mesh with run output alone, rewrite the other
# meshes' configs, and refresh the run-tree tools (scripts replaced by rename, so a running chain keeps its own copy).
set -e
SRC=$(dirname $(readlink -f $0))
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
OLD=SST_p0j8.015_p0inf3.247; NEW=SST_p0j7.4945_p0inf3.247
if [ "$UPDATE" = 1 ]; then
  MESHES=""
  for m in Medium Fine FineJet FineJet2 VeryFine; do
    if [ -n "$(find $R/${NEW}_mesh$m -name history.csv -not -path '*/test/*' 2>/dev/null)" ]; then echo "$m has run output: not touching it"
    else MESHES="$MESHES $m"; fi
  done
  [ -n "$MESHES" ] || { echo "nothing to update"; exit 1; }
  for f in $SRC/*.sh $SRC/README.md; do cp $f $R/${NEW}_tools/.new_$(basename $f); mv $R/${NEW}_tools/.new_$(basename $f) $R/${NEW}_tools/$(basename $f); done
  echo "updated from cfd-tools $(git -C $SRC rev-parse HEAD) on $(date) (UPDATE=1)" >> $R/${NEW}_tools/COMMIT
else
  MESHES="Coarse Medium Fine FineJet FineJet2 VeryFine"
  for m in $MESHES; do
    [ -e $R/${NEW}_mesh$m ] && { echo "$R/${NEW}_mesh$m exists: not touching it"; exit 1; }
  done
  [ -e $R/${NEW}_tools ] && { echo "$R/${NEW}_tools exists"; exit 1; }
  # tools: copies from the clone (watcher and plots from the SST study), with the commit they came from
  mkdir $R/${NEW}_tools
  cp $SRC/../${OLD}_tools/{watch_fj.sh,box_snap.py,box_change.py,plot_conv_roi.py,stage_summary.py} $SRC/*.sh $SRC/README.md $R/${NEW}_tools/
  echo "copied from cfd-tools $(git -C $SRC rev-parse HEAD) on $(date)" > $R/${NEW}_tools/COMMIT
fi
declare -A FROM=([Medium]=Coarse [Fine]=Medium [FineJet]=Fine [FineJet2]=FineJet [VeryFine]=FineJet2)

INLET='MARKER_INLET= ( inlet_jet, 288.15, 749446.0, 1.0, 0.0, 0.0 )   % p0j 7.4945 bar, T0j 288.15 K: from p = 7.0588 bar, V = 99.11 m/s (user, 2026-11-02)'
NOTE='% Phase 2: inlet_jet is a total-conditions inlet, p0j = 7.4945 bar, T0j = 288.15 K, flow in +x (corrected, user 2026-11-02).'
edit() {  # edit <header line> <source cfg> <mesh file>: the common changes, to stdout
  sed -e "1s|.*|$1|" \
      -e "s|$OLD|$NEW|g" \
      -e "s|^% Nothing stops on convergence. Made by.*|% Nothing stops on convergence. Made by ../../../${NEW}_tools/setup_p0j7.4945.sh from the $OLD configs.|" \
      -e "s|^% Phase 2: inlet_jet is.*|$NOTE|" \
      -e "s|^MARKER_INLET=.*|$INLET|" \
      -e "s|^MESH_FILENAME=.*|MESH_FILENAME= $3|" $2
}
cd $R
for m in $MESHES; do
  if [ $m = Coarse ]; then tmpl=${OLD}_meshMedium; last=t4_cap3; else tmpl=${OLD}_mesh$m; last=t3_cap3; fi
  mesh=$(sed -n 's/^MESH_FILENAME= *//p' ${OLD}_mesh$m/2_jet_staged/$last/$last.cfg)
  if [ $m = Coarse ]; then
    [ -s ${OLD}_meshCoarse/2_jet_staged/t4_cap3/restart.dat ] || { echo "missing ${OLD}_meshCoarse/2_jet_staged/t4_cap3/restart.dat"; exit 1; }
    prev="from $OLD meshCoarse t4_cap3 (8.015 bar inlet), same mesh"
  else
    g=${NEW}_mesh$m/0_interp/interp.cfg; mkdir -p $(dirname $g)
    edit "% srp_M2.5, case $NEW, mesh$m, k-omega SST (user, 2026-11-03). Interpolation of mesh${FROM[$m]} t3_cap3 (1 iteration at CFL 1e-8)" \
      ${OLD}_mesh$m/0_interp/interp.cfg $mesh > $g
    sed -i "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $R/${NEW}_mesh${FROM[$m]}/2_jet_staged/t3_cap3/restart|" $g
    prev="from mesh${FROM[$m]} t3_cap3, interpolated (0_interp)"
  fi
  for s in t1_cap0.5 t2_cap1 t3_cap3; do
    g=${NEW}_mesh$m/2_jet_staged/$s/$s.cfg; mkdir -p $(dirname $g)
    edit "% srp_M2.5, case $NEW, mesh$m, k-omega SST (user, 2026-11-02). Stage $s: $prev" $tmpl/2_jet_staged/$s/$s.cfg $mesh > $g
    if [ $s = t1_cap0.5 ]; then
      [ $m = Coarse ] && src=$R/${OLD}_meshCoarse/2_jet_staged/t4_cap3/restart || src=../../0_interp/restart
      sed -i "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $src|" $g
    fi
    prev="from $s"
  done
done
echo "built:"; find $(for m in $MESHES; do echo ${NEW}_mesh$m; done) -name "*.cfg" -not -path "*/test/*" | sort
echo "inlet lines:"; grep -h "^MARKER_INLET=\|^INLET_TYPE" $(find ${NEW}_mesh* -name "*.cfg" -not -path "*/test/*") | sort | uniq -c
echo "mesh files:"; for m in $MESHES; do echo "  $m $(grep -h '^MESH_FILENAME' $(find ${NEW}_mesh$m -name '*.cfg' -not -path '*/test/*') | sort -u)"; done
echo "sources:"; grep -H "^SOLUTION_FILENAME" $(find ${NEW}_mesh* -name "*.cfg" -not -path "*/test/*" | sort)
echo "old-name references left:"; grep -n "$OLD" $(find ${NEW}_mesh* -name "*.cfg" -not -path "*/test/*") | grep -v "^[^:]*:[13]:\|SOLUTION_FILENAME" || echo none
