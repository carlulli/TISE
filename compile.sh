IDIR="./include"
MDIR="./modules"

gcc -o main.o -I ${IDIR} ${MDIR}/geometry.c ${MDIR}/hamiltonian.c ${MDIR}/powermethod.c ${MDIR}/observables.c ${MDIR}/linearalgebra.c ${MDIR}/conjugategradient.c ./main.c
