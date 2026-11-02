#!/bin/bash
# Build SST_p0j7.4945_p0inf3.247_* (user, 2026-11-02): the SST mesh study (SST_p0j8.015_p0inf3.247_*) repeated with the
# corrected jet inlet. The user gave static p = 7.0588 bar and V = 99.11 m/s at the jet inlet (plenum face); with
# T0 = 288.15 K (option A): T = 288.15 - 99.11^2/(2*1004.703) = 283.26 K, M = 0.2937, p0 = 7.0588*(1 + 0.2 M^2)^3.5
# = 7.4945 bar. The inlet stays TOTAL_CONDITIONS; only the p0j value changes.
# Each mesh restarts from its own SST 8.015 bar final solution (coarse t4_cap3, the others t3_cap3; same mesh, no
# interpolation), then t1 cap 0.5 (5000) / t2 cap 1 (25000) / t3 cap 3 (20000), as the SST meshes after coarse. The coarse
# configs are the meshMedium ones with the coarse mesh file. Nothing stops on convergence.
# Run on M3 from the cfd-tools clone: bash <clone>/srp_M2.5/runs/SST_p0j7.4945_p0inf3.247_tools/setup_p0j7.4945.sh
set -e
SRC=$(dirname $(readlink -f $0))
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
OLD=SST_p0j8.015_p0inf3.247; NEW=SST_p0j7.4945_p0inf3.247
for m in Coarse Medium Fine FineJet FineJet2 VeryFine; do
  [ -e $R/${NEW}_mesh$m ] && { echo "$R/${NEW}_mesh$m exists: not touching it"; exit 1; }
done
[ -e $R/${NEW}_tools ] && { echo "$R/${NEW}_tools exists"; exit 1; }

# tools: copies from the clone (watcher and plots from the SST study), with the commit they came from
mkdir $R/${NEW}_tools
cp $SRC/../${OLD}_tools/{watch_fj.sh,box_snap.py,box_change.py,plot_conv_roi.py,stage_summary.py} $SRC/*.sh $SRC/README.md $R/${NEW}_tools/
echo "copied from cfd-tools $(git -C $SRC rev-parse HEAD) on $(date)" > $R/${NEW}_tools/COMMIT

INLET='MARKER_INLET= ( inlet_jet, 288.15, 749446.0, 1.0, 0.0, 0.0 )   % p0j 7.4945 bar, T0j 288.15 K: from p = 7.0588 bar, V = 99.11 m/s (user, 2026-11-02)'
NOTE='% Phase 2: inlet_jet is a total-conditions inlet, p0j = 7.4945 bar, T0j = 288.15 K, flow in +x (corrected, user 2026-11-02).'
cd $R
for m in Coarse Medium Fine FineJet FineJet2 VeryFine; do
  if [ $m = Coarse ]; then tmpl=${OLD}_meshMedium; last=t4_cap3; else tmpl=${OLD}_mesh$m; last=t3_cap3; fi
  mesh=$(sed -n 's/^MESH_FILENAME= *//p' ${OLD}_mesh$m/2_jet_staged/$last/$last.cfg)
  [ -s ${OLD}_mesh$m/2_jet_staged/$last/restart.dat ] || { echo "missing ${OLD}_mesh$m/2_jet_staged/$last/restart.dat"; exit 1; }
  prev="from $OLD mesh$m $last (8.015 bar inlet), same mesh"
  for s in t1_cap0.5 t2_cap1 t3_cap3; do
    g=${NEW}_mesh$m/2_jet_staged/$s/$s.cfg; mkdir -p $(dirname $g)
    sed -e "1s|.*|% srp_M2.5, case $NEW, mesh$m, k-omega SST (user, 2026-11-02). Stage $s: $prev|" \
        -e "s|$OLD|$NEW|g" \
        -e "s|^% Nothing stops on convergence. Made by.*|% Nothing stops on convergence. Made by ../../../${NEW}_tools/setup_p0j7.4945.sh from the $OLD configs.|" \
        -e "s|^% Phase 2: inlet_jet is.*|$NOTE|" \
        -e "s|^MARKER_INLET=.*|$INLET|" \
        -e "s|^MESH_FILENAME=.*|MESH_FILENAME= $mesh|" \
        $tmpl/2_jet_staged/$s/$s.cfg > $g
    [ $s = t1_cap0.5 ] && sed -i "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $R/${OLD}_mesh$m/2_jet_staged/$last/restart|" $g
    prev="from $s"
  done
done
echo "built:"; find ${NEW}_mesh* -name "*.cfg" | sort
echo "inlet lines:"; grep -h "^MARKER_INLET=\|^INLET_TYPE" $(find ${NEW}_mesh* -name "*.cfg") | sort | uniq -c
echo "mesh files:"; for m in Coarse Medium Fine FineJet FineJet2 VeryFine; do echo "  $m $(grep -h '^MESH_FILENAME' ${NEW}_mesh$m/2_jet_staged/*/*.cfg | sort -u)"; done
echo "t1 sources:"; grep -H "^SOLUTION_FILENAME" ${NEW}_mesh*/2_jet_staged/t1_cap0.5/t1_cap0.5.cfg
echo "old-name references left:"; grep -n "$OLD" $(find ${NEW}_mesh* -name "*.cfg") | grep -v "^[^:]*:1:\|SOLUTION_FILENAME" || echo none
