#!/bin/bash
# Build SSTsarkar_p0j8.015_p0inf3.247_* (user, 2026-10-31): the SST run (SST_p0j8.015_p0inf3.247_*) repeated with the
# Sarkar compressibility correction. Every config is a copy of the SST one as it was run (including the Phase 1 CFL of
# start 0.1, cap 0.5, floor 0.1); only SST_OPTIONS and the run-folder names change.
set -e
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
OLD=SST_p0j8.015_p0inf3.247; NEW=SSTsarkar_p0j8.015_p0inf3.247
for m in Coarse Medium Fine FineJet FineJet2 VeryFine; do
  [ -e $R/${NEW}_mesh$m ] && { echo "$R/${NEW}_mesh$m exists: not touching it"; exit 1; }
done
[ -e $R/${NEW}_tools ] && { echo "$R/${NEW}_tools exists"; exit 1; }
mkdir $R/${NEW}_tools
cp $R/${OLD}_tools/{watch_fj.sh,box_snap.py,box_change.py,plot_conv_roi.py,stage_summary.py,setup_sst.sh} $R/${NEW}_tools/
cd $R
for f in $(find ${OLD}_mesh* -maxdepth 3 -name "*.cfg" -not -path "*failed*" -not -path "*test*" | sort); do
  g=$NEW${f#$OLD}; mkdir -p $(dirname $g)
  { echo "% SSTsarkar (user, 2026-10-31): copy of ../$OLD/${f#*/} with the Sarkar compressibility correction; nothing else changed."
    sed -e "s|$OLD|$NEW|g" \
        -e 's|^SST_OPTIONS=.*|SST_OPTIONS= V2003m, COMPRESSIBILITY-SARKAR   % SST-2003m with the Sarkar compressibility correction (user, 2026-10-31)|' $f
  } > $g
done
echo "built:"; find ${NEW}_mesh* -name "*.cfg" | sort
echo "SST_OPTIONS lines:"; grep -h "^SST_OPTIONS" $(find ${NEW}_mesh* -name "*.cfg") | sort | uniq -c
echo "old-name references left:"; grep -l "$OLD" $(find ${NEW}_mesh* -name "*.cfg") || echo none
echo "interp sources:"; grep -h "^SOLUTION_FILENAME" ${NEW}_mesh*/0_interp/interp.cfg
