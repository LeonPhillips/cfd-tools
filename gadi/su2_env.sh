# Source this to use SU2 v8.5.0 on Gadi:  source /scratch/os85/lp3540/su2/SU2-v8.5.0/su2_env.sh
# (copied there by su2_setup_login.sh). Same build options as the M3 install.
module load gcc/13.2.0 intel-mpi/2021.13.1 python3/3.11.7 >/dev/null 2>&1
SU2_BASE=/scratch/os85/lp3540/su2
source $SU2_BASE/su2-venv/bin/activate
export SU2_HOME=$SU2_BASE/SU2-v8.5.0/src
export SU2_RUN=$SU2_BASE/SU2-v8.5.0/install/bin
export PATH=$SU2_RUN:$PATH
export PYTHONPATH=$SU2_RUN:${PYTHONPATH:-}
