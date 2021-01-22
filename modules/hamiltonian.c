#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "geometry.h"
#include "linearalgebra.h"

/*parameters of the Hamiltonian that need to be set*/
static int N;
static double mass;
static double *V;
static double min = 0.0; //minimum value of the potential;
static int shift = 0;


/* Sets the mass and initializes the V array in accordance to the N passed by the user*/
void set_kinetic_params(double m) {
	mass = m;
	static int temp = 0;
	if(temp==0 ) {
	N = get_N();
	V = malloc(sizeof(double)*N);
	if(V == NULL)
    {
        printf("Error! memory not allocated.");
    }

	temp++;
}
	else{
		printf("Hamiltonian already computed\n");
	}

}
/* unfortunate attempt to avoid using set_potential  */
/*void set_potential(const char* fname, union Potential potential ){
	if strcmp(fname, "ZERO") == 0 {
		set_zero_potential();
	} else if strcmp(fname, "HARMONIC") == 0{
		set_harmonic_potential(potential.k);
	} else if strcmp(fname, "WELL") == 0{
		set_well_potential(potential.eps, potential.a, potential.b);
	} else if strcmp(fname, "WALL") == 0{
		set_wall_potential(potential.eps, potential.a, potential.b);
	} else if strcmp(fname, "COULOMB") == 0{
		set_coulomb_potential();
	} else {
		fprintf( stderr,"Error unknown option '%s', expected one of ZERO, HARMONIC, WELL, WALL, COULOMB.", fname);
	}
}
*/

/* sets the potential tied to a number passed by the user */

void set_coulomb_potential() {
	int i;
	for (i = 1 ; i < N / 2 ; i++) {
		V[N / 2 + i] = - 1 / (i);
		V[N / 2 - i] = - 1 / (i);
	}
	V[N / 2] = V[N / 2 + 1];
}


void set_zero_potential() {
	int i;
	for(i = 0; i < N; i++) {
		V[i] = 0;
	}
}


void set_harmonic_potential() {
	double	k = 1;
	int i;
	for(i = 0 ; i < N ; i++)	{

	V[i] = 0.5*k * (i - (N-1)/2)*(i-(N-1)/2);
	}
}

void set_well_potential() {
	double eps = 1;
	int a = (int)N/4;
	for(int j = 0; j < N; j++) {
		V[j] = eps;
		while(j > a && j < N-a) {
			V[j] = 0;
			j +=1;
		}
	}
}

void set_wall_potential() {
	double eps = 1;
	for(int j = 0; j < N; j++) {
		V[j] = 0;
		while(j > (N+1)/2 && j < N) {
			V[j] =  eps;
			j += 1;
		}
	}
}
/* makes the hamiltonian definite positive by shifting V by the lowest point */

void set_minV() {
	double minV = 0;
	for(int i = 0; i < N; i++) {
		if(V[i] < minV) {
			minV = V[i];
		}
	}
	min = minV;
}

void set_potential(int pot) {

	if(pot == 0) {
		set_zero_potential();
		printf("[hamiltonian.c | set_potential()] Chosen POTENTIAL: Zero-potential \n");
	}
	else if(pot == 1) {
		set_harmonic_potential();
		printf("[hamiltonian.c | set_potential()] Chosen POTENTIAL: harmonic potential \n");
	}
	else if(pot == 2 ) {
		set_well_potential();
		printf("[hamiltonian.c | set_potential()] Chosen POTENTIAL: well potential \n");
	}
	else if(pot == 3) {
		set_wall_potential();
		printf("[hamiltonian.c | set_potential()] Chosen POTENTIAL: wall potential \n");
	}
	else {
		printf("[hamiltonian.c | set_potentilal()] No potential assigned! Possibel choices 0,1,2,3\n");
		exit(-1);
	}
	set_minV();
}

void set_Hdefpos() {
	if (shift == 0) {
		for(int j = 0; j < N; j++) {
		V[j] -= min;
	}
	shift++;
	printf("The Hamiltonian has been made definite positive \n");
}
	else {
		printf("The hamiltonian is already definite positive \n");
	}
}

double get_minV() {
	if (min == 0) {
		printf("[hamiltonian.c | get_minV()] min is still 0 (to check min=%f)\n", min);
	}
	return min;
}

void H(double complex *in, double complex *out) {

	double complex delta[N];
	int i;

	for (i = 1; i < (N - 1); i++) {//finite difference
		delta[i] = in[i - 1] - in[i] + in[i + 1] - in[i];

	}
	delta[0] = in[1] -2*in[0]; //Boundary conditions
	delta[N-1] = in[N-2] -2*in[N-1];


	for (i = 0; i < N; i++) {//applies H to input vector

		 out[i] = - 1 / (2 * mass) * (delta[i]) + V[i] * in[i]; //writes into the output vector

	 }
}

void H_defpos(double complex *in, double complex *out) {
	if(shift == 0) {
	set_Hdefpos();
	}
	H(in,out);
}



double average_state_energy(double complex *psi) {
	double complex Hpsi[N];
	H(psi,Hpsi);

	return creal(scalar_product(psi, Hpsi,N) / norm(psi,N)*norm(psi,N));


}

double  average_kinetic_energy(double complex *psi) {
	double complex delta[N];
	int i;

	for (i = 1; i < (N - 1); i++) {//finite difference
		delta[i] = psi[i - 1] - psi[i] + psi[i + 1] - psi[i];

	}
	delta[0] = psi[1] -2*psi[0]; //Boundary conditions
	delta[N-1] = psi[N-2] -2*psi[N-1];

 for( i = 0; i < N; i++) {

	 delta[i] = 1 / ( 2 * mass) * delta[i];
  }

	return creal(scalar_product(psi,delta,N) / norm(psi,N)*norm(psi,N));

	}

	double average_potential_energy(double complex *psi) {
		double complex Vpsi[N];
		for( int i = 0; i < N; i ++) {
			Vpsi[i] = V[i]*psi[i];
		}
		return creal(scalar_product(psi,Vpsi,N)/norm(psi,N)*norm(psi,N));
	}

void print_hamiltonian_info() {
	printf("Hamiltonian parameters: \nmass = %f\nN = %d\n", mass, N);
	for(int i = 0; i < N; i++){
	printf("Potential V[%d] = %f\n",i,V[i]);
	}
}

void free_V() { free(V);}
