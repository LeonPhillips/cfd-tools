#!/bin/bash
# Gadi login-node part of the SU2 v8.5.0 install (needs internet, which compute nodes don't have):
# clone the source, make the Python venv, and run meson configure (it fetches subprojects).
# The compile itself is build_su2.pbs. Run: bash su2_setup_login.sh   (from the cfd-tools clone or piped over ssh)
set -e
B=/scratch/os85/lp3540/su2
mkdir -p $B/SU2-v8.5.0
module load gcc/13.2.0 intel-mpi/2021.13.1 python3/3.11.7

# Python venv (numpy/matplotlib for the analysis tools, mpi4py against this Intel MPI for pysu2)
if [ ! -x $B/su2-venv/bin/python ]; then
  python3 -m venv --system-site-packages $B/su2-venv
  $B/su2-venv/bin/pip install -q --upgrade pip
  MPICC=mpicc $B/su2-venv/bin/pip install -q --no-binary mpi4py mpi4py
  $B/su2-venv/bin/pip install -q numpy scipy matplotlib
fi
$B/su2-venv/bin/python -c "import numpy, matplotlib, mpi4py; print('venv ok', numpy.__version__, mpi4py.__version__)"

# Source
[ -d $B/SU2-v8.5.0/src/.git ] || git clone -q --depth 1 --branch v8.5.0 https://github.com/su2code/SU2.git $B/SU2-v8.5.0/src
cd $B/SU2-v8.5.0/src
git describe --tags

# Configure (same options as M3)
source $B/su2-venv/bin/activate
CC=mpicc CXX=mpicxx ./meson.py setup build --prefix=$B/SU2-v8.5.0/install --buildtype=release \
  -Dcustom-mpi=true -Dwith-mpi=enabled -Dwith-omp=true -Denable-pywrapper=true
cp "$(dirname "$(readlink -f "$0")")/su2_env.sh" $B/SU2-v8.5.0/su2_env.sh 2>/dev/null || true
echo "configure done $(date)"
