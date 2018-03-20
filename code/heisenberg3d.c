#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

/*
Requirements: GNU Scientific Library (gsl)

gcc -lgsl heisenberg3d.c
*/

#define ROWS 10
#define COLS 10

typedef struct {
  double x;
  double y;
  double z;
} spin_t;

typedef struct {
  spin_t layer1[ROWS][COLS];
  spin_t layer2[ROWS][COLS];
} lattice_t;

gsl_rng * rng;


void initialize_lattice(lattice_t*);
void gen_random_spin(spin_t*);
void simulate(lattice_t*, int, double);
void sweep(lattice_t*, double);



int main(){

  int rows, cols;
  double J_intra, J_inter, B, k1, k2, init_T;
  lattice_t lattice;

  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (rng, 0);


  initialize_lattice(&lattice);

}

void initialize_lattice(lattice_t* lattice){
  int i, j, k;
  for(i = 0; i < ROWS; i++){
    for(j = 0; j < COLS; j++){
      gen_random_spin(&(*lattice).layer1[i][j]);
      gen_random_spin(&(*lattice).layer2[i][j]);
    }
  }
}

void gen_random_spin(spin_t* spin){
    double x1, x2, mag_sq;

    x1 = gsl_rng_uniform(rng)*2.0 - 1.0;
    x2 = gsl_rng_uniform(rng)*2.0 - 1.0;
    mag_sq = gsl_pow_2(x1) + gsl_pow_2(x2);

    while(mag_sq >=1){
      x1 = gsl_rng_uniform(rng)*2.0 - 1.0;
      x2 = gsl_rng_uniform(rng)*2.0 - 1.0;
      mag_sq = gsl_pow_2(x1) + gsl_pow_2(x2);
    }

    spin->x = 2.0*x1*sqrt(1.0-mag_sq);
    spin->y = 2.0*x2*sqrt(1.0-mag_sq);
    spin->z = 1.0-2.0*mag_sq;

    //printf("x: %f\n", spin->x);
    //printf("y: %f\n", spin->y);
    //printf("z: %f\n", spin->z);

}

void simulate(lattice_t* lattice, int num_sweeps, double T){
  int i;
  for(i = 0; i < num_sweeps; i++){
    sweep(lattice, T);
  }
}

void sweep(lattice_t* lattice, double T){


}
