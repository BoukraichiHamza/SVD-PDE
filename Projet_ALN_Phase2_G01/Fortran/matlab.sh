#! /bin/sh
#
# Scrip t de paramétrage d'environnement pour Matlab
#
make clean
make
make matlabsetup
export PATH=/applications/matlab/bin:$PATH
make mexfile
export OMP_NUM_THREADS=4
export LD_PRELOAD=/usr/lib/gcc/x86_64-linux-gnu/4.9/libgomp.so
matlab
