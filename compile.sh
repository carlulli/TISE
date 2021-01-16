################################################################################
# This script can be used to compile the main.c when inside TISE directory
# It creates an executable main.o that needs to be called with the parameter
# for the TISE simulation
# like this:
# ./main.o [int NUM] [double mass] [int potential number] [double tol] [double res]
################################################################################

IDIR="./include"
MDIR="./modules"

gcc -o main.o -I ${IDIR} ${MDIR}/geometry.c ${MDIR}/hamiltonian.c ${MDIR}/powermethod.c ${MDIR}/observables.c ${MDIR}/linearalgebra.c ${MDIR}/conjugategradient.c ./main.c
