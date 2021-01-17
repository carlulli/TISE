#!/usr/bin/bash
IDIR="../include"
MDIR="../modules"
gcc -Wall -o anal -I $IDIR $MDIR/geometry.c $MDIR/assert.c $MDIR/hamiltonian.c $MDIR/linearalgebra.c analytic.c -lm
