#!/usr/bin/bash
IDIR="../include"
MDIR="../modules"
gcc -Wall -o lin -I $IDIR $MDIR/geometry.c $MDIR/assert.c $MDIR/hamiltonian.c $MDIR/linearalgebra.c $MDIR/wavefunction.c linearity.c -lm
echo "type ./lin [N] [MASS] [SET POTENTIAL] [tolerance]"
