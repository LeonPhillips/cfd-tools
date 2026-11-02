#!/bin/bash
# Interpolate the meshFineNozzle n3_cap3 solution onto meshVeryFine with SU2_CFD's own restart interpolation.
# Run on the desktop node (not the login node): bash <this folder>/interp.sh   (24 ranks, an even number)
source /scratch/qk75/LeonPhillips/PhD/cfd/su2/SU2-v8.5.0/su2_env.sh
cd /scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/runs/SA_p0j8.015_p0inf3.247_meshVeryFine/0_interp || exit 1
{ echo "host $(hostname)  ranks $(nproc)"; date; } > interp.log
mpirun -n $(nproc) SU2_CFD -t 1 interp.cfg >> interp.log 2>&1
echo "SU2_CFD exit=$?  $(date)" >> interp.log
