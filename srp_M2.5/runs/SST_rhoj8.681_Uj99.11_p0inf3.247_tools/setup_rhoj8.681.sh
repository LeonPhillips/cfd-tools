#!/bin/bash
# Build SST_rhoj8.681_Uj99.11_p0inf3.247_mesh{Coarse,FineJet} (user, 2026-11-02): SST with a MASS_FLOW jet inlet.
# The user gave static p = 7.0588 bar and V = 99.11 m/s at the jet inlet (plenum face, r 0-2.8787 mm, x -8.124 mm).
# With T0 = 288.15 K: T = 283.26 K, rho = 705880/(287.058*283.26) = 8.681 kg/m3, mass flow ~0.0224 kg/s (full 360 deg).
# SU2 MASS_FLOW reads MARKER_INLET= ( marker, density, velocity magnitude, direction ) (CEulerSolver::BC_Inlet).
# meshCoarse from scratch with the long staging of the original SST study: Phase 1 (nozzle plugged, 20000, CFL 0.1 -> cap
# 0.5, floor 0.1), then t1 cap 0.3 / t2 cap 0.5 / t3 cap 1 / t4 cap 3, 50000 each. meshFineJet: coarse t4 interpolated
# onto it (0_interp), then t1 cap 0.5 (5000) / t2 cap 1 (25000) / t3 cap 3 (20000). Configs are copies of the
# SST_p0j8.015_p0inf3.247 ones as run; only the inlet, the names and the interpolation source change.
# Run on M3 from the cfd-tools clone: bash <clone>/srp_M2.5/runs/SST_rhoj8.681_Uj99.11_p0inf3.247_tools/setup_rhoj8.681.sh
set -e
SRC=$(dirname $(readlink -f $0))
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
OLD=SST_p0j8.015_p0inf3.247; NEW=SST_rhoj8.681_Uj99.11_p0inf3.247
for m in Coarse FineJet; do
  [ -e $R/${NEW}_mesh$m ] && { echo "$R/${NEW}_mesh$m exists: not touching it"; exit 1; }
done
[ -e $R/${NEW}_tools ] && { echo "$R/${NEW}_tools exists"; exit 1; }

mkdir $R/${NEW}_tools
cp $SRC/../${OLD}_tools/{watch_fj.sh,box_snap.py,box_change.py,plot_conv_roi.py,stage_summary.py} $SRC/*.sh $SRC/README.md $R/${NEW}_tools/
echo "copied from cfd-tools $(git -C $SRC rev-parse HEAD) on $(date)" > $R/${NEW}_tools/COMMIT

INLET='MARKER_INLET= ( inlet_jet, 8.681, 99.11, 1.0, 0.0, 0.0 )   % MASS_FLOW: rho 8.681 kg/m3, V 99.11 m/s (p 7.0588 bar, T 283.26 K; user, 2026-11-02)'
NOTE='% Phase 2: inlet_jet is a MASS_FLOW inlet, rho = 8.681 kg/m3, V = 99.11 m/s, flow in +x (user, 2026-11-02).'
cd $R
for f in $(find ${OLD}_meshCoarse ${OLD}_meshFineJet -maxdepth 3 -name "*.cfg" -not -path "*failed*" -not -path "*test*" | sort); do
  g=$NEW${f#$OLD}; mkdir -p $(dirname $g)
  { echo "% $NEW (user, 2026-11-02): copy of ../$OLD/${f#*/} with the MASS_FLOW jet inlet. Made by ${NEW}_tools/setup_rhoj8.681.sh."
    sed -e "s|$OLD|$NEW|g" \
        -e "s|^INLET_TYPE=.*|INLET_TYPE= MASS_FLOW|" \
        -e "s|^% Phase 2: inlet_jet is.*|$NOTE|" \
        -e "s|^MARKER_INLET=.*|$INLET|" $f
  } > $g
done
# meshFineJet starts from this case's coarse t4, not from fine
sed -i "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $R/${NEW}_meshCoarse/2_jet_staged/t4_cap3/restart|" ${NEW}_meshFineJet/0_interp/interp.cfg
sed -i "1s|.*|% srp_M2.5, case $NEW, meshFineJet, k-omega SST (user, 2026-11-02). Stage t1_cap0.5: from meshCoarse t4_cap3, interpolated (0_interp)|" \
  ${NEW}_meshFineJet/2_jet_staged/t1_cap0.5/t1_cap0.5.cfg
echo "built:"; find ${NEW}_mesh* -name "*.cfg" | sort
echo "inlet lines:"; grep -h "^MARKER_INLET=\|^INLET_TYPE" $(find ${NEW}_mesh* -name "*.cfg") | sort | uniq -c
echo "phase 1 jet:"; grep -h "inlet_jet" ${NEW}_meshCoarse/1_freestream/phase1.cfg
echo "sources:"; grep -H "^SOLUTION_FILENAME" $(find ${NEW}_mesh* -name "*.cfg" | sort)
echo "old-name references left:"; grep -n "$OLD" $(find ${NEW}_mesh* -name "*.cfg") | grep -v "^[^:]*:1:" || echo none
