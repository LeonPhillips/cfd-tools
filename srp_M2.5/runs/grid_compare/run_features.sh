#!/bin/bash
# Key-feature positions for every final srp_M2.5 solution (user, 2026-11-04): SST p0j 7.4945 bar (5 meshes),
# SST MASS_FLOW rho 8.681 (coarse, fineJet), and the p0j 8.015 bar SST, SST-Sarkar and SA mesh studies.
# Final stage = the run's only 2_jet_staged/*_cap3. Resamples a 0.025 mm box with Mach, Density, Velocity, then runs
# feature_positions.py. Run on a compute node (e.g. the desktop job): bash run_features.sh [label-filter regex]
# Output: $R/grid_compare/features/ (box_features_all_0.025mm.npz, feature_positions_all.csv, features.log)
# MODE=machdisk TAG=sst bash run_features.sh '^SST': 0.005 mm box x 4-5, r 0-0.5 mm and machdisk_fine.py -> machdisk_fine_sst.csv
TOOLS=$(dirname $(readlink -f $0))
R=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs
O=$R/grid_compare/features; mkdir -p $O; cd $O
FILTER=${1:-.}
exec > >(tee -a features.log) 2>&1
echo "== start $(date) on $(hostname), cfd-tools $(git -C $TOOLS rev-parse --short HEAD), filter '$FILTER'"
ARGS=()
add() {  # add <run folder prefix> <model> <inlet> <meshes...>
  local p=$1 mo=$2 in=$3; shift 3
  for m in "$@"; do
    local f=( $R/${p}_mesh$m/2_jet_staged/*_cap3/flow.vtu )
    [ ${#f[@]} -eq 1 ] && [ -s ${f[0]} ] || { echo "SKIP $p mesh$m: ${#f[@]} *_cap3/flow.vtu"; continue; }
    local st=$(basename $(dirname ${f[0]})) lab=${mo}__${in}__${m}__
    lab=$lab$st
    [[ $lab =~ $FILTER ]] && ARGS+=( "$lab=${f[0]}" )
  done
}
add SST_p0j7.4945_p0inf3.247 SST p0j7.4945 Coarse Medium Fine FineJet FineJet2
add SST_rhoj8.681_Uj99.11_p0inf3.247 SST massflow_rhoj8.681 Coarse FineJet
add SST_p0j8.015_p0inf3.247 SST p0j8.015 Coarse Medium Fine FineJet FineJet2 VeryFine
add SSTsarkar_p0j8.015_p0inf3.247 SSTsarkar p0j8.015 Coarse Medium Fine FineJet FineJet2 VeryFine
add SA_p0j8.015_p0inf3.247 SA p0j8.015 Coarse Medium Fine FineNozzle FineJet FineJet2 VeryFine
echo "${#ARGS[@]} solutions"; printf '  %s\n' "${ARGS[@]}"
TAG=${TAG:-$([ "$FILTER" = . ] && echo all || echo test)}
if [ "$MODE" = machdisk ]; then   # MODE=machdisk: fine box around the Mach disk only (machdisk_fine.py)
  NPZ=box_machdisk_${TAG}_0.005mm.npz; export BOX=4.0,5.0,0.0,0.5 H_MM=0.005 VARS=Mach,Density
else
  NPZ=box_features_${TAG}_0.025mm.npz; export VARS=Mach,Density,Velocity H_MM=0.025
fi
( unset LD_PRELOAD; module load paraview/5.13 >/dev/null 2>&1
  pvbatch --mesa $TOOLS/resample_box.py $NPZ "${ARGS[@]}" 2>&1 | grep -vi warning )
[ -s $NPZ ] || { echo "RESAMPLE FAILED"; exit 1; }
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh >/dev/null
if [ "$MODE" = machdisk ]; then
  python $TOOLS/machdisk_fine.py $NPZ machdisk_fine_${TAG}.csv
else
  python $TOOLS/feature_positions.py $NPZ feature_positions_${TAG}.csv
fi
ls -l $O
echo "== done $(date)"
