# TISE
This code calculates the grounds state and its energy of the time independent Schrödinger equation in 1D.

Additionally, the following observables are calculated: average position and momentum, standard deviation of position and momentum.

The code is mainly written in c, with shell scripts for compiling.

After downloading the code, the local directory structure has to be as follows:
Inside TISE:
  - include directory: contains all necessary header files
  - modules directory: contains the necessary c code files of each module
  - test directory: contains the test files that test each module
  - output directory: the directory where the simulation output files are written
  - main.c: The main.c file to compile prior to running simulations.
  - compile.sh: a shell script that compiles the main.c code

### Running the code
After downloading and saving the code structure, run ``` compile.sh ``` from within the TISE directory to compile the code.
(Reminder: include compiling code for the test codes.)

Following successful compiling run ``` main.o [NUM] [mass] [potential number] [tol] [res] ```
The parameters are as follows:
 - NUM = Number of lattice points for the discretized volume. This needs to be an integer > 0 and odd.
 - mass = mass of the discretized hamiltonian. This will be a double.
 - potential number = 0 zero potential, 1 harmonic potential, 2 well potential, 3 wall potential
 - tol = tolerance for the power method. Will be a double.
 - res = residue for the conjugate gradient. Will be a double.

 Running the testing: ??

### bout the modules

- asser.c:

- conjugategradient.c: Method that calculates M^-1x out of initial Mx=b with input matrix M and vector b

- geometry.c: sets the geometry parameters, here sets N = NUM

- hamiltonian.c: sets kin part of hamiltonian with mass and N and the potential, contains all potentails, makes hamiltonina positive definite to use in power method

- linearalgebra.c: contains functions that calculate scalar product and norm

- observables.c: contains functions that calculate obseravbles average position and momentum, standard deviation of position and momentum.

- powermethod.c: a method that calculates the largest eigenvalue of given matrix. Takes as input a matrix (in this case the inverse hamiltonian) and output vector. Writes eigenstate as output vector and returns eigenvalue.

### About the testing

- analytic.c:
- geometrytest.c:
- hermitian.c:
- power_algo.c:
- test_conjugate.c:
