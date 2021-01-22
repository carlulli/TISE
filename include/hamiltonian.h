#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <complex.h>


/*sets the mass parameter that is used in the hamiltonian*/
void set_kinetic_params(double m);

//void set_potential(const char* fname, union Potential potential);

/*******************************************************************************
The following functions set the potential used by the hamiltonian.c module. The
following options are available

1 set_zero_potential()              : V(n) = 0
2 set_harmonic_potential(double k)  : V(n) = k/2 * (n-(N-1)/2)^2
3 set_well_potential(eps,a,b)       :
4 set_wall_potential(eps,a,b)       :
5 set_coulomb_potential();
The potential must be set exactly once by the main program.
*******************************************************************************/
void set_potential(int pot);
void set_zero_potential();
void set_harmonic_potential();
void set_well_potential();
void set_wall_potential();
void set_coulomb_potential();

void set_minV();
void set_Hdefpos();
double get_minV();
/*******************************************************************************
This function applies the Hamiltonian to the wave function in and returns the
result in the wavefunction out.

in, out = pointers to arrays with N elements, where N is the number of points of
          the lattice; the memory for these arrays is assumed to be already
          allocated
*******************************************************************************/
void H_defpos(double complex *in, double complex *out);
void H(double complex *in, double complex *out);

double average_state_energy(double complex *psi);
double average_kinetic_energy(double complex *psi);
double average_potential_energy(double complex *psi);
/*******************************************************************************
This function prints the parameters of the Hamiltonian on standard output.
*******************************************************************************/
void print_hamiltonian_info();

/*free the allocated memory of V*/
void free_V();

#endif /* HAMILTONIAN_H */
