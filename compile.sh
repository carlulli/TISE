################################################################################
# This script can be used to compile the main.c when inside TISE directory
# It creates an executable main.o that needs to be called with the parameter
# It also compiles all test files and creates their executables.
# The test executables are created inside the test directory.
################################################################################

IDIR="./include"
MDIR="./modules"
TDIR="./test"
SDIR="./scripts"

gcc -o ${SDIR}/main.o -I ${IDIR} ${MDIR}/geometry.c ${MDIR}/hamiltonian.c ${MDIR}/powermethod.c ${MDIR}/observables.c ${MDIR}/linearalgebra.c ${MDIR}/conjugategradient.c main.c -lm
gcc -o ${TDIR}/analytic.o -I ${IDIR} ${MDIR}/geometry.c ${MDIR}/hamiltonian.c ${MDIR}/linearalgebra.c ${MDIR}/assert.c ${TDIR}/analytic.c -lm
gcc -o ${TDIR}/geometrytest.o -I ${IDIR} ${MDIR}/geometry.c ${TDIR}/geometrytest.c -lm
gcc -o ${TDIR}/hermitian.o -I ${IDIR} ${MDIR}/geometry.c ${MDIR}/hamiltonian.c ${MDIR}/wavefunction.c ${MDIR}/linearalgebra.c ${MDIR}/assert.c ${TDIR}/hermitian.c -lm
gcc -o ${TDIR}/power_algo.o -I ${IDIR} ${MDIR}/geometry.c ${MDIR}/powermethod.c ${TDIR}/power_algo.c -lm
gcc -o ${TDIR}/test_conjugate.o ${TDIR}/test_conjugate.c -lm
