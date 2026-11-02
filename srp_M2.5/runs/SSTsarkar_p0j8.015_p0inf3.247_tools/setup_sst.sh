#!/bin/bash
# Build the SST (k-omega SST) run folders for srp_M2.5, case p0j 8.015 bar / p0inf 3.247 bar (user, 2026-10-30).
# Configs are the SA ones with only the turbulence model, turbulence inflow and the run-control lines changed.
# Usage (login or desktop node; it only writes text files): bash setup_sst.sh
#   SST_options: V2003m, no compressibility correction; freestream TI 1 %, mu_t/mu 1; jet inlet TI 2 %, mu_t/mu 10.
#   coarse: Phase 1 (freestream, nozzle plugged, fixed 20000 iterations, CFL 1 -> 50), then jet on:
#           t1 cap 0.3 (first order 500 iterations, second order by 1000), t2 cap 0.5, t3 cap 1, t4 cap 3; 50000 each.
#   others: interpolated from the previous mesh's last stage, then t1 cap 0.5 (5000), t2 cap 1 (25000), t3 cap 3 (20000).
set -e
B=/scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5; R=$B/runs; M=$B/mesh
TOOLS=$R/SST_p0j8.015_p0inf3.247_tools
P1=$R/SA_p0j8.015_p0inf3.247_meshCoarse/1_freestream/phase1.cfg            # SA Phase 1
T=$R/SA_p0j8.015_p0inf3.247_meshFineJet2/2_jet_staged/j3_cap3/j3_cap3.cfg   # SA stage (probes, snapshots, no-stop)
I=$R/SA_p0j8.015_p0inf3.247_meshFineJet2/0_interp/interp.cfg                # SA interpolation (1 iteration, CFL 1e-8)
run() { echo $R/SST_p0j8.015_p0inf3.247_mesh$1; }
declare -A MESH=([Coarse]=srp_M2.5_coarse_v4 [Medium]=srp_M2.5_meshMedium [Fine]=srp_M2.5_meshFine
                 [FineJet]=srp_M2.5_meshFineJet [FineJet2]=srp_M2.5_meshFineJet2 [VeryFine]=srp_M2.5_meshVeryFine)
# mesh <- source mesh for the interpolation
declare -A FROM=([Medium]=Coarse [Fine]=Medium [FineJet]=Fine [FineJet2]=FineJet [VeryFine]=Fine)
for m in Coarse Medium Fine FineJet FineJet2 VeryFine; do [ -e $(run $m) ] && { echo "$(run $m) exists: not touching it"; exit 1; }; done
for m in "${MESH[@]}"; do [ -s $M/$m.su2 ] || { echo "missing mesh $m"; exit 1; }; done

# SA -> SST physics (all configs); jet-on configs also get the jet-inlet turbulence
SST=(-e 's|^KIND_TURB_MODEL= SA$|KIND_TURB_MODEL= SST|' -e '/^% standard SA-neg/d'
     -e 's|^SA_OPTIONS=.*|SST_OPTIONS= V2003m   % standard SST-2003m, no compressibility correction (user, 2026-10-30)|'
     -e 's|^FREESTREAM_NU_FACTOR=.*|FREESTREAM_TURBULENCEINTENSITY= 0.01   % freestream TI 1 % (user, 2026-10-30)\nFREESTREAM_TURB2LAMVISCRATIO= 1.0      % freestream mu_t/mu (user, 2026-10-30)|'
     -e 's|RMS_NU_TILDE|RMS_TKE, RMS_DISSIPATION|')
JET=(-e 's|^\(MARKER_INLET= .*\)$|\1\nMARKER_INLET_TURBULENT= ( inlet_jet, 0.02, 10.0 )   % jet TI 2 %, mu_t/mu 10 (user, 2026-10-30)|')
NOHEAD=(-e '1,/^% ---- physics/{/^% ---- physics/!d}')

# stage <run> <name> <solution> <CFL start> <cap> <iterations> <ramp: YES|NO> <text>
stage() {
  local d=$1/2_jet_staged/$2; mkdir -p $d
  { echo "% srp_M2.5, case SST_p0j8.015_p0inf3.247, mesh${1##*_mesh}, k-omega SST (user, 2026-10-30). Stage $2: $8"
    echo "% CFL starts at $4, ramps x1.01 per iteration to a cap of $5, halved on trouble, floor 0.1; $6 iterations, full length."
    echo "% Nothing stops on convergence. Made by ../../../SST_p0j8.015_p0inf3.247_tools/setup_sst.sh from the SA configs."
    sed "${NOHEAD[@]}" "${SST[@]}" "${JET[@]}" \
        -e "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $3|" -e "s|^CFL_NUMBER=.*|CFL_NUMBER= $4|" \
        -e "s|^CFL_ADAPT_PARAM=.*|CFL_ADAPT_PARAM= ( 0.5, 1.01, 0.1, $5, 0.001, 0 )|" -e "s|^ITER=.*|ITER= $6|" \
        -e "s|^MESH_FILENAME=.*|MESH_FILENAME= $M/${MESH[${1##*_mesh}]}.su2|" \
        $T
  } > $d/$2.cfg
  if [ "$7" = YES ]; then   # jet switched on: first order for 500 iterations, second order by 1000 (as SA s1)
    sed -i -e 's|^% fully second order:.*|% first order for 500 iterations, then a linear ramp to fully second order at iteration 1000|' \
           -e 's|^RAMP_MUSCL= NO$|RAMP_MUSCL= YES\nRAMP_MUSCL_COEFF= ( 500, 10, 500 )\nKIND_MUSCL_RAMP= ITERATION|' $d/$2.cfg
  fi
}
# interp <run> <source restart>
interp() {
  local d=$1/0_interp; mkdir -p $d
  { echo "% Interpolation of $2 onto mesh${1##*_mesh} (k-omega SST, 2026-10-30): 1 iteration at CFL 1e-8."
    echo "% SU2_CFD interpolates the restart itself because its point count differs from the mesh (needs an even number of ranks)."
    sed "${NOHEAD[@]}" "${SST[@]}" "${JET[@]}" -e "s|^SOLUTION_FILENAME=.*|SOLUTION_FILENAME= $2|" \
        -e "s|^MESH_FILENAME=.*|MESH_FILENAME= $M/${MESH[${1##*_mesh}]}.su2|" $I
  } > $d/interp.cfg
}

# coarse: Phase 1 + four 50000-iteration stages
C=$(run Coarse); mkdir -p $C/1_freestream
{ echo "% srp_M2.5, case SST_p0j8.015_p0inf3.247, mesh v4 (meshCoarse), k-omega SST (user, 2026-10-30). PHASE 1: freestream M2.5,"
  echo "% jet off (nozzle plugged), from a uniform freestream start. As the SA Phase 1 except the turbulence model and the"
  echo "% length: a fixed 20000 iterations (nothing stops on convergence). First order for 500 iterations, then second order"
  echo "% by iteration 1500; CFL 1 -> 50."
  sed "${NOHEAD[@]}" "${SST[@]}" -e "s|^CONV_FIELD=.*|CONV_FIELD= RMS_DENSITY|" \
      -e "s|^CONV_RESIDUAL_MINVAL=.*|CONV_RESIDUAL_MINVAL= -20   % unreachable: runs never stop on a convergence test (user, 2026-10-28)|" \
      -e "s|^ITER=.*|ITER= 20000|" -e "s|^MESH_FILENAME=.*|MESH_FILENAME= $M/${MESH[Coarse]}.su2|" $P1
} > $C/1_freestream/phase1.cfg
stage $C t1_cap0.3 ../../1_freestream/restart 0.1 0.3 50000 YES "jet switched on, from the Phase 1 restart"
stage $C t2_cap0.5 ../t1_cap0.3/restart 0.3 0.5 50000 NO "from t1_cap0.3"
stage $C t3_cap1   ../t2_cap0.5/restart 0.5 1.0 50000 NO "from t2_cap0.5"
stage $C t4_cap3   ../t3_cap1/restart   1.0 3.0 50000 NO "from t3_cap1"
# the other meshes
for m in Medium Fine FineJet FineJet2 VeryFine; do
  d=$(run $m); src=$(run ${FROM[$m]})
  [ ${FROM[$m]} = Coarse ] && last=t4_cap3 || last=t3_cap3
  interp $d $src/2_jet_staged/$last/restart
  stage $d t1_cap0.5 ../../0_interp/restart 0.1 0.5 5000 NO "from mesh${FROM[$m]} $last, interpolated (0_interp)"
  stage $d t2_cap1 ../t1_cap0.5/restart 0.5 1.0 25000 NO "from t1_cap0.5"
  stage $d t3_cap3 ../t2_cap1/restart 1.0 3.0 20000 NO "from t2_cap1"
done
echo "built:"; for m in Coarse Medium Fine FineJet FineJet2 VeryFine; do find $(run $m) -name "*.cfg" | sort; done
