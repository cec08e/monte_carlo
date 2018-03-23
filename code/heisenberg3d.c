#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>

/*
Requirements: GNU Scientific Library (gsl) and CBLAS

gcc -lgsl -lgslcblas heisenberg3d.c
*/

#define ROWS 10       /* Number of rows in each lattice layer */
#define COLS 10       /* Number of columns in each lattice layer */
#define RADIUS .4     /* Radius of tangent disc in perturbing function */
#define K1 -1
#define K2 -1
#define J_INTRA 1
#define J_INTER 1
#define B_EXT 1

typedef struct {
  double x;
  double y;
  double z;
} spin_t;             /* Spin type structure - x, y, and z components of spin */

/*
typedef struct {
  spin_t layer1[ROWS][COLS];
  spin_t layer2[ROWS][COLS];
} lattice_t;
*/

typedef spin_t lattice_t[2][ROWS][COLS];

gsl_rng * rng;


void initialize_lattice(lattice_t);
void gen_random_spin(spin_t*);
void simulate(lattice_t, int, double);
void sweep(lattice_t, double);
void perturb_spin(spin_t* , spin_t*);
double calc_delta_E(lattice_t, spin_t* , spin_t*, int, int, int);




int main(){

  int rows, cols;
  double J_intra, J_inter, B, k1, k2, init_T;
  lattice_t lattice;

  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (rng, 0);


  initialize_lattice(lattice);
  simulate(lattice, 1, .1);

}

void initialize_lattice(lattice_t lattice){
  int i, j, k;
  for(i = 0; i < ROWS; i++){
    for(j = 0; j < COLS; j++){
      gen_random_spin(&lattice[0][i][j]);
      gen_random_spin(&lattice[1][i][j]);
    }
  }
}

void gen_random_spin(spin_t* spin){
    double x1, x2, mag_sq;

    x1 = gsl_rng_uniform(rng)*2.0 - 1.0;
    x2 = gsl_rng_uniform(rng)*2.0 - 1.0;
    mag_sq = gsl_pow_2(x1) + gsl_pow_2(x2);

    while(mag_sq >= 1){
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

void simulate(lattice_t lattice, int num_sweeps, double T){
  int i;
  for(i = 0; i < num_sweeps; i++){
    sweep(lattice, T);
  }
}

void sweep(lattice_t lattice, double T){
  int i,num;
  int layer, row, col;
  double delta_E;
  double random_num;
  spin_t temp_spin;
  /* Perform as many MC steps as there are lattice sites */
  for(i=0; i < 2*ROWS*COLS; i++){
    /* Choose a random spin site on a random lattice */
    layer = gsl_rng_get(rng) % 2;
    row = gsl_rng_get(rng) % ROWS;
    col = gsl_rng_get(rng) % COLS;
    printf("Chosen site on layer %d: (%d,%d)\n", layer, row, col);
    /*
    if(layer == 0 && row == 0 && col == 7){
      scanf("%d", &num);

    }
    */

    /* Generate a perturbation of the spin at the chosen lattice site */
    perturb_spin(&temp_spin, &lattice[layer][row][col]);
    delta_E = calc_delta_E(lattice, &temp_spin, &lattice[layer][row][col], layer, row, col);
    random_num = gsl_rng_uniform(rng);
    printf("random number: %f\n", random_num);
    if( !( (delta_E > 0) && (random_num >= gsl_sf_exp(-(1.0/T)*delta_E)) ) ){
      lattice[layer][row][col].x = temp_spin.x;
      lattice[layer][row][col].y = temp_spin.y;
      lattice[layer][row][col].z = temp_spin.z;
      printf("Accepted. \n");
    }


  }

}

void perturb_spin(spin_t* temp_spin, spin_t* spin){
  double r, arg;
  double x, y, z, x_new, y_new, u_x, u_y;
  double theta, phi;
  double norm;
  gsl_vector* spin_vector = gsl_vector_alloc(3);
  gsl_vector_set(spin_vector, 0, spin->x);
  gsl_vector_set(spin_vector, 1, spin->y);
  gsl_vector_set(spin_vector, 2, spin->z);

  gsl_vector* point_vector = gsl_vector_alloc(3);
  gsl_vector* temp_vector = gsl_vector_calloc(3);

  gsl_matrix* rot_matrix = gsl_matrix_alloc(3,3);

  /* Generate random position on the tangent disc */
  r = sqrt(gsl_rng_uniform(rng))*RADIUS;
  arg = gsl_rng_uniform(rng)*2*M_PI;

  /* Express random point as cartesian coordinate in same reference frame as spin */
  x = r*sin(arg);
  y = r*cos(arg);
  z = 0;

  /* Grab spherical coordinates of spin vector */

  theta = acos(spin->z);
  if(spin->x == 0){
    if(spin->y > 0)
      phi = M_PI/2.0;
    else
      phi = 3*M_PI/2.0;
  }
  else
    phi = atan(spin->y/spin->x);

  /* Rotate random point with phi */
  x_new = x*cos(phi) - y*sin(phi);
  y_new = x*sin(phi) + y*cos(phi);

  gsl_vector_set(point_vector, 0, x_new);
  gsl_vector_set(point_vector, 1, y_new);
  gsl_vector_set(point_vector, 2, z);

  /* Now, rotate random point with theta - rotate around y' = -sin(phi), cos(phi) axis */
  u_x = -sin(phi);
  u_y = cos(phi);

  /* Creating rotation matrix */
  gsl_matrix_set(rot_matrix, 0, 0, cos(theta) +  gsl_pow_2(u_x)*(1-cos(theta)));
  gsl_matrix_set(rot_matrix, 0, 1, u_x*u_y*(1-cos(theta)));
  gsl_matrix_set(rot_matrix, 0, 2, u_y*sin(theta));
  gsl_matrix_set(rot_matrix, 1, 0, u_x*u_y*(1-cos(theta)));
  gsl_matrix_set(rot_matrix, 1, 1, cos(theta) + gsl_pow_2(u_y)*(1-cos(theta)));
  gsl_matrix_set(rot_matrix, 1, 2, -u_x*sin(theta));
  gsl_matrix_set(rot_matrix, 2, 0, -u_y*sin(theta));
  gsl_matrix_set(rot_matrix, 2, 1, u_x*sin(theta));
  gsl_matrix_set(rot_matrix, 2, 2,  cos(theta));


  gsl_blas_dgemv(CblasNoTrans,1.0, rot_matrix, point_vector,0.0,temp_vector);
  gsl_vector_add(temp_vector, spin_vector);
  norm = sqrt(gsl_pow_2(gsl_vector_get(temp_vector, 0))
              + gsl_pow_2(gsl_vector_get(temp_vector, 1))
              + gsl_pow_2(gsl_vector_get(temp_vector, 2)));
  gsl_vector_scale(temp_vector, 1/norm);

  temp_spin->x = gsl_vector_get(temp_vector, 0);
  temp_spin->y = gsl_vector_get(temp_vector, 1);
  temp_spin->z = gsl_vector_get(temp_vector, 2);

  printf("Generated perturbation: %f, %f, %f \n", temp_spin->x, temp_spin->y, temp_spin->z);

  gsl_vector_free(spin_vector);
  gsl_vector_free(point_vector);
  gsl_vector_free(temp_vector);
  gsl_matrix_free(rot_matrix);
}

double calc_delta_E(lattice_t lattice, spin_t* temp_spin, spin_t* spin, int layer, int row, int col){

  double delta_dot_neighbor;
  double delta_dot_inter;
  double delta_a;
  double delta_E;
  int i;

  /* Calculate change in spin */
  gsl_vector* delta_vector = gsl_vector_alloc(3);
  gsl_vector_set(delta_vector, 0, temp_spin->x - spin->x);
  gsl_vector_set(delta_vector, 1, temp_spin->y - spin->y);
  gsl_vector_set(delta_vector, 2, temp_spin->z - spin->z);

  /* Calculate neighbor sum */
  gsl_vector* neighbor_vector = gsl_vector_alloc(3);
  gsl_vector_set(neighbor_vector, 0, lattice[layer][(((row-1)%ROWS) + ROWS) % ROWS][col].x
                                    + lattice[layer][(((row+1)%ROWS) + ROWS) % ROWS][col].x
                                    + lattice[layer][row][(((col-1)%COLS) + COLS) % COLS].x
                                    + lattice[layer][row][(((col+1)%COLS) + COLS) % COLS].x);
  gsl_vector_set(neighbor_vector, 1, lattice[layer][(((row-1)%ROWS) + ROWS) % ROWS][col].y
                                    + lattice[layer][(((row+1)%ROWS) + ROWS) % ROWS][col].y
                                    + lattice[layer][row][(((col-1)%COLS) + COLS) % COLS].y
                                    + lattice[layer][row][(((col+1)%COLS) + COLS) % COLS].y);
  gsl_vector_set(neighbor_vector, 2, lattice[layer][(((row-1)%ROWS) + ROWS) % ROWS][col].z
                                    + lattice[layer][(((row+1)%ROWS) + ROWS) % ROWS][col].z
                                    + lattice[layer][row][(((col-1)%COLS) + COLS) % COLS].z
                                    + lattice[layer][row][(((col+1)%COLS) + COLS) % COLS].z);

  //printf("Neighbor z: \n");
  //printf("row-1 mod rows: %d \n", (((row-1)%ROWS) + ROWS) % ROWS);
  //printf("%f, %f, %f, %f \n", lattice[layer][(((row-1)%ROWS) + ROWS) % ROWS][col].z, lattice[layer][(((row+1)%ROWS) + ROWS) % ROWS][col].z, lattice[layer][row][(col-1)%COLS].z, lattice[layer][row][(col+1)%COLS].z);
  gsl_blas_ddot(delta_vector, neighbor_vector, &delta_dot_neighbor);


  gsl_vector* inter_vector = gsl_vector_alloc(3);
  gsl_vector_set(inter_vector, 0, lattice[(layer+1)%2][row][col].x);
  gsl_vector_set(inter_vector, 1, lattice[(layer+1)%2][row][col].y);
  gsl_vector_set(inter_vector, 2, lattice[(layer+1)%2][row][col].z);

  gsl_blas_ddot(delta_vector, inter_vector, &delta_dot_inter);


  /* Calculate anisotropy change */
  if(layer == 0)
    delta_a = K1*(gsl_pow_2(gsl_vector_get(delta_vector,2)) + 2*gsl_vector_get(delta_vector,2)*spin->z);
  else
    delta_a = K2*(gsl_pow_2(gsl_vector_get(delta_vector,2)) + 2*gsl_vector_get(delta_vector,2)*spin->z);

  /*
  for (i = 0; i < 3; i++){
    printf ("neighbor_vector_%d = %g\n", i, gsl_vector_get (neighbor_vector, i));
  }
  printf("delta_dot_neighbor: %f\n", delta_dot_neighbor);
  printf("delta_dot_inter: %f\n", delta_dot_inter);
  printf("delta_a: %f\n", delta_a);
  */

  delta_E = -J_INTRA*delta_dot_neighbor + J_INTER*delta_dot_inter + delta_a - B_EXT*gsl_vector_get(delta_vector,2);
  printf("Delta E is %f \n", delta_E);

  return delta_E;



}
