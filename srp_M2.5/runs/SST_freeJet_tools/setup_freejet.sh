#!/bin/bash
# Build the free-jet run folders on Gadi (user, 2026-11-23) under /scratch/os85/lp3540/su2/cases/srp_M2.5_freeJet/runs:
#   SST_p0j6.6_pa101325_freeJet_mesh{Coarse,Medium,Fine}/
#     0_interp/interp.cfg          (Medium, Fine) the previous mesh's t4_cap3 restart, interpolated by SU2_CFD
#     0_seed/                      (Coarse) base.txt -> SEED_BASE (default $C/seed/base_meshCoarse_restart.dat, a
#                                  coarse restart with the nozzle running); run_fj_chain.sh writes restart.dat from it
#                                  with seed_jet.py: a decaying round jet along the axis out to the outlet (user,
#                                  2026-11-23, option 1)
#     2_jet_staged/t1_cap0.5_c001  first 1000-iteration chunk of t1, CFL 0.1 -> cap 0.5 (coarse: from 0_seed,
#                                  FIRST ORDER (MUSCL_FLOW= NO) in every t1 chunk (user, 2026-11-23) to damp
#                                  the starting vortex; coarse t2 ramps to second order by its iteration 800. Medium,
#                                  fine: from 0_interp, second order)
#     t1_chunk.cfg                 template for the next t1 chunks (CFL 0.5, cap 0.5, from the chunk before;
#                                  run_fj_chain.sh fills PREV_CHUNK). t1 runs chunk by chunk until the jet reaches the
#                                  outlet (axis Ux > 20 m/s at x = 220 mm), then one more chunk (user, 2026-11-23).
#                                  t1_final -> the last chunk. Iterations per mesh: t1 (n x 1000) + 49,000.
#                  t2_cap1    1000, CFL 0.5 -> cap 1, from t1_final
#                  t3_cap3    8000, CFL 1 -> cap 3
#                  t4_cap3   40000, CFL 1 -> cap 3          (50,000 iterations per mesh)
#     title.txt
# Usage (from the Gadi clone): bash setup_freejet.sh           Refuses to touch a run folder that exists.
set -e
T=$(dirname $(readlink -f $0))
C=${C:-/scratch/os85/lp3540/su2/cases/srp_M2.5_freeJet}
CASE=SST_p0j6.6_pa101325_freeJet
COMMIT=$(git -C $T rev-parse --short HEAD)
declare -A FROM=([Medium]=Coarse [Fine]=Medium)
for m in Coarse Medium Fine; do
  [ -s $C/mesh/srp_M2.5_freeJet_mesh$m.su2 ] || { echo "missing mesh $C/mesh/srp_M2.5_freeJet_mesh$m.su2"; exit 1; }
  [ -e $C/runs/${CASE}_mesh$m ] && { echo "$C/runs/${CASE}_mesh$m exists: not touching it"; exit 1; }
done
mk() {  # mk <out cfg> <mesh> <stage note> <restart> <solution> <muscl ramp> <cfl0> <adapt> <cap> <iter> <wrt> <overwrite>
  mkdir -p $(dirname $1)
  sed -e "s|@CASE@|$CASE|" -e "s|@MESH@|$2|" -e "s|@STAGE_NOTE@|$3|" -e "s|@COMMIT@|$COMMIT|" \
      -e "s|@RESTART@|$4|" -e "s|@SOLUTION@|$5|" -e "s|@MUSCL_RAMP@|$6|" -e "s|@CFL0@|$7|" -e "s|@CFL_ADAPT@|$8|" \
      -e "s|@CAP@|$9|" -e "s|@ITER@|${10}|" -e "s|@WRT@|${11}|" -e "s|@OVERWRITE@|${12}|" \
      -e "s|@MUSCL@|${MUSCL:-YES}|" -e "s|@MESHFILE@|$C/mesh/srp_M2.5_freeJet_mesh$2.su2|" $T/freejet_template.cfg > $1
  ! grep -n "@[A-Z_0-9]*@" $1 || { echo "unfilled placeholder in $1"; exit 1; }
}
NORAMP='RAMP_MUSCL= NO   % fully second order'
RAMP='RAMP_MUSCL= YES   % cold start: first order for 300 iterations, then a linear ramp to second order by iteration 800
RAMP_MUSCL_COEFF= ( 300, 10, 500 )
KIND_MUSCL_RAMP= ITERATION'
RAMP=$(printf '%s' "$RAMP" | sed ':a;N;$!ba;s/\n/\\n/g')
RAMP2='RAMP_MUSCL= YES   % coarse t2: linear ramp from the first-order t1 to second order by iteration 800
RAMP_MUSCL_COEFF= ( 0, 10, 800 )
KIND_MUSCL_RAMP= ITERATION'
RAMP2=$(printf '%s' "$RAMP2" | sed ':a;N;$!ba;s/\n/\\n/g')
for m in Coarse Medium Fine; do
  R=$C/runs/${CASE}_mesh$m
  if [ $m = Coarse ]; then
    SEED_BASE=${SEED_BASE:-$C/seed/base_meshCoarse_restart.dat}
    [ -s $SEED_BASE ] || { echo "missing seed base $SEED_BASE"; exit 1; }
    mkdir -p $R/0_seed; echo $SEED_BASE > $R/0_seed/base.txt
    sol=../../0_seed/restart; rs=YES; ramp=$NORAMP; t2ramp=$RAMP2; mu1=NO
    from="from the seeded jet column (0_seed, seed_jet.py), first order"
  else
    mk $R/0_interp/interp.cfg $m "Interpolation of mesh${FROM[$m]} t4_cap3 (1 iteration at CFL 1e-8; SU2_CFD interpolates the restart itself because its point count differs from the mesh, even number of ranks)" \
       YES $C/runs/${CASE}_mesh${FROM[$m]}/2_jet_staged/t4_cap3/restart "$NORAMP" 1e-8 NO 3.0 1 1 YES
    sol=../../0_interp/restart; rs=YES; ramp=$NORAMP; t2ramp=$NORAMP; mu1=YES; from="from mesh${FROM[$m]} t4_cap3, interpolated (0_interp)"
  fi
  MUSCL=$mu1 mk $R/2_jet_staged/t1_cap0.5_c001/t1_cap0.5_c001.cfg $m "Stage t1_cap0.5, chunk 1: 1000 iterations $from; CFL starts at 0.1, ramps x1.01 to a cap of 0.5, halved on trouble, floor 0.1." \
     $rs $sol "$ramp" 0.1 YES 0.5 1000 "1000, 1000, 1000" YES
  MUSCL=$mu1 mk $R/t1_chunk.cfg $m "Stage t1_cap0.5, next chunk: 1000 iterations from PREV_CHUNK; CFL 0.5, cap 0.5, halved on trouble, floor 0.1; MUSCL_FLOW $mu1." \
     YES ../PREV_CHUNK/restart "$NORAMP" 0.5 YES 0.5 1000 "1000, 1000, 1000" YES
  mk $R/2_jet_staged/t2_cap1/t2_cap1.cfg $m "Stage t2_cap1: 1000 iterations from the last t1 chunk (t1_final); CFL 0.5 -> cap 1." \
     YES ../../t1_final/restart "$t2ramp" 0.5 YES 1.0 1000 "1000, 1000, 1000" NO
  mk $R/2_jet_staged/t3_cap3/t3_cap3.cfg $m "Stage t3_cap3: 8000 iterations from t2_cap1; CFL 1 -> cap 3." \
     YES ../t2_cap1/restart "$NORAMP" 1.0 YES 3.0 8000 "1000, 1000, 1000" NO
  mk $R/2_jet_staged/t4_cap3/t4_cap3.cfg $m "Stage t4_cap3: 40000 iterations from t3_cap3; CFL 1 -> cap 3." \
     YES ../t3_cap3/restart "$NORAMP" 1.0 YES 3.0 40000 "1000, 1000, 1000" NO
  echo "SST free jet, p0j 6.6 bar into 1 atm, mesh$m, SLAU2-VW" > $R/title.txt
  echo "made by SST_freeJet_tools/setup_freejet.sh, cfd-tools $COMMIT, $(date)" > $R/SOURCE.txt
done
cd $C/runs
echo "built:"; find ${CASE}_mesh* -name "*.cfg" | sort
echo "check:"; for f in $(find ${CASE}_mesh* -name "*.cfg" | sort); do
  printf '%-60s %s | %s | %s | %s\n' $f "$(grep -h '^RESTART_SOL' $f)" "$(grep -h '^SOLUTION_FILENAME' $f | sed 's|.*/runs/||')" \
    "$(grep -h '^ITER' $f)" "$(grep -h '^CFL_ADAPT_PARAM' $f | cut -d, -f4)"; done
