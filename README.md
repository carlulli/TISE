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
After downloading and saving the code structure, run ``` build.sh ``` from within the TISE directory to compile the code.
This compiles and links the main code and the testing files. 

Following successful compiling run ``` main.o [NUM] [mass] [potential number] [tol] [res] ```
The parameters are as follows:
 - NUM = Number of lattice points for the discretized volume. This needs to be an integer > 0 and odd.
 - mass = mass of the discretized hamiltonian. This will be a double.
 - potential number = 0 zero potential, 1 harmonic potential, 2 well potential, 3 wall potential
 - tol = tolerance for the power method. Will be a double.
 - res = residue for the conjugate gradient. Will be a double.

 Running the testing: After running ``` build.sh ``` the test output files should be created inside the test directory. From there all tests can be run.

### bout the modules

- assert.c: Check that the number of parameter input in the bash is the correct one

- conjugategradient.c: Method that calculates M^-1x out of initial Mx=b with input matrix M and vector b

- geometry.c: sets the geometry parameters, here sets N = NUM

- hamiltonian.c: sets kin part of hamiltonian with mass and N and the potential, contains all potentails, makes hamiltonina positive definite to use in power method

- linearalgebra.c: contains functions that calculate scalar product and norm

- observables.c: contains functions that calculate obseravbles average position and momentum, standard deviation of position and momentum.

- powermethod.c: a method that calculates the largest eigenvalue of given matrix. Takes as input a matrix (in this case the inverse hamiltonian) and output vector. Writes eigenstate as output vector and returns eigenvalue.

### About the testing

- analytic.c: Calculates the eigenspectrum for the free particle case with the hamiltonian.c module
- geometrytest.c:tests that the geometrical parameters are correctly stored
- hermitian.c: tests hermitianity of the hamiltonian.c module
- power_algo.c: tests the power algorithm for a general def pos matrix
- test_conjugate.c: tests the conjugate algorthm for a general def pos matrix
