#!/bin/bash
# Build the meshFineJet2 run folder from the meshVeryFine configs (2026-10-29). Run on the desktop node.
# Usage: bash setup_fineJet.sh <run folder> <source restart, without .dat>   (relative to <run folder>/0_interp)
# Stage configs are made from the very-fine v3 config by make_cfg (only the lines below change), plus:
#   probes (CUSTOM_OUTPUTS, Mach at the points in PROBES) in the history, and numbered flow_<iter>.vtu snapshots
#   every SNAP iterations (WRT_VOLUME_OVERWRITE= NO) for the box-change check.
B=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5
V=$B/runs/SA_p0j8.015_p0inf3.247_meshVeryFine
FJ=${1:?run folder}; SRC=${2:?source restart}
SNAP=${SNAP:-1000}
T=$V/2_jet_staged/v3_cap3/v3_cap3.cfg
grep -q "^WRT_VOLUME_OVERWRITE\|^CUSTOM_OUTPUTS" $T && { echo "template already sets these: fix make_cfg"; exit 1; }
# Mach = |velocity| / sound speed (probes cannot use MACH directly); x [m], r [m]: see plot_conv_roi.py for what each one is
M="sqrt(pow(VELOCITY_X, 2) + pow(VELOCITY_Y, 2)) / SOUND_SPEED"
PROBES="p1 : Probe{$M}[0.00400, 0.00000]; p2 : Probe{$M}[0.00750, 0.00000]; p3 : Probe{$M}[0.00828, 0.00070];\
 p4 : Probe{$M}[0.00752, 0.00148]; p5 : Probe{$M}[0.00724, 0.00200]; p6 : Probe{$M}[0.00680, 0.00230];\
 p7 : Probe{$M}[0.00524, 0.00300]"
mkdir -p $FJ/0_interp $FJ/2_jet_staged

# make_cfg <out> <solution> <CFL start> <cap> <iterations> <header text>
make_cfg() {
  { echo "% srp_M2.5, case SA_p0j8.015_p0inf3.247, meshFineJet2, PHASE 2 STAGED CFL (user, 2026-10-29). $6"
    echo "% CFL starts at $3, ramps x1.01 per iteration to a cap of $4, halved on trouble, floor 0.1; up to $5 iterations."
    echo "% Each stage runs its full ITER (user, 2026-10-29: no early exit). Nothing stops on convergence."
    sed -e '1,/^$/d' \
        -e "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $2|" \
        -e "s|^CFL_NUMBER=.*|CFL_NUMBER= $3|" \
        -e "s|^CFL_ADAPT_PARAM=.*|CFL_ADAPT_PARAM= ( 0.5, 1.01, 0.1, $4, 0.001, 0 )|" \
        -e "s|^ITER=.*|ITER= $5|" \
        -e "s|^MESH_FILENAME=.*|MESH_FILENAME= $B/mesh/srp_M2.5_meshFineJet2.su2|" \
        -e "s|^\(HISTORY_OUTPUT=.*\) )|\1, CUSTOM )|" \
        -e "s|^OUTPUT_WRT_FREQ=.*|OUTPUT_WRT_FREQ= 1000, $SNAP, 1000   % restart, flow_<iter>.vtu snapshot, surface|" $T
    echo
    echo "% ---- region-of-interest monitoring (2026-10-29) ----"
    echo "WRT_VOLUME_OVERWRITE= NO"
    echo "CUSTOM_OUTPUTS= '$PROBES'"
  } > $1
}
# stage name, CFL start, cap, iterations (as meshFineJet j1-j3)
STAGES="j1_cap0.5 0.1 0.5 5000
j2_cap1 0.5 1.0 25000
j3_cap3 1.0 3.0 20000"
prev=../../0_interp/restart
while read s c0 cap n; do
  mkdir -p $FJ/2_jet_staged/$s
  make_cfg $FJ/2_jet_staged/$s/$s.cfg $prev $c0 $cap $n "Stage $s, from ${prev#../}."
  prev=../$s/restart
done <<< "$STAGES"

# interpolation: the very-fine interp.cfg with the source and mesh swapped (1 iteration at CFL 1e-8)
sed -e '1,2d' -e "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $SRC|" \
    -e "s|^MESH_FILENAME=.*|MESH_FILENAME= $B/mesh/srp_M2.5_meshFineJet2.su2|" $V/0_interp/interp.cfg |
  { echo "% Interpolation of $SRC onto meshFineJet2 (2026-10-29): 1 iteration at CFL 1e-8."
    echo "% SU2_CFD interpolates the restart itself because its point count differs from the mesh (needs an even number of ranks)."
    cat; } > $FJ/0_interp/interp.cfg
cat > $FJ/0_interp/interp.sh <<EOF
#!/bin/bash
# Interpolate $SRC onto meshFineJet2 with SU2_CFD's own restart interpolation.
# Run on the desktop node (not the login node): bash <this folder>/interp.sh   (24 ranks, an even number)
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
cd $FJ/0_interp || exit 1
{ echo "host \$(hostname)  ranks \$(nproc)"; date; } > interp.log
mpirun -n \$(nproc) SU2_CFD -t 1 interp.cfg >> interp.log 2>&1
echo "SU2_CFD exit=\$?  \$(date)" >> interp.log
EOF
ls -R $FJ | head -20
