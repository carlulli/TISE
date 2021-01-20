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

gcc -Wall -o main.o -I ${IDIR} ${MDIR}/geometry.c ${MDIR}/hamiltonian.c ${MDIR}/powermethod.c ${MDIR}/observables.c ${MDIR}/linearalgebra.c ${MDIR}/conjugategradient.c main.c -lm
gcc -Wall -o ${TDIR}/analytic.o -I ${IDIR} ${MDIR}/geometry.c ${MDIR}/hamiltonian.c ${MDIR}/linearalgebra.c ${MDIR}/assert.c ${TDIR}/analytic.c -lm
gcc -Wall -o ${TDIR}/geometrytest.o -I ${IDIR} ${MDIR}/geometry.c ${TDIR}/geometrytest.c -lm
gcc -Wall -o ${TDIR}/hermitian.o -I ${IDIR} ${MDIR}/geometry.c ${MDIR}/hamiltonian.c ${MDIR}/wavefunction.c ${MDIR}/linearalgebra.c ${MDIR}/assert.c ${TDIR}/hermitian.c -lm
gcc -Wall -o ${TDIR}/power_algo.o -I ${IDIR} ${MDIR}/geometry.c ${MDIR}/powermethod.c ${MDIR}/assert.c  ${MDIR}/linearalgebra.c  ${TDIR}/power_algo.c -lm
gcc -Wall -o ${TDIR}/test_conjugate.o -I ${IDIR} ${MDIR}/geometry.c ${MDIR}/assert.c ${MDIR}/linearalgebra.c ${MDIR}/conjugategradient.c ${TDIR}/test_conjugate.c -lm
