#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <time.h>

/*
Requirements: GNU Scientific Library (gsl) and CBLAS

gcc -lgsl -lgslcblas heisenberg3d.c
ALT: gcc -fPIC -shared -o heisenberg3d.so -lgsl -lgslcblas heisenberg3d.c

sudo gdb python3
run plot_h3d.py
*/

#define ROWS 10       /* Number of rows in each lattice layer */
#define COLS 10       /* Number of columns in each lattice layer */
#define RADIUS .4   /* Radius of tangent disc in perturbing function */
#define J_INTRA 1
#define J_INTER 1
#define INIT_T 5

double B_EXT = 1;
double K1 = -1;
double K2 = -1;

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
lattice_t lattice;


void initialize_lattice();
void gen_random_spin(spin_t*);
void simulate( int, double);
int sweep( double);
void perturb_spin(spin_t* , spin_t*);
double calc_delta_E( spin_t* , spin_t*, int, int, int);
double calc_magnetization( int);
void M_v_B();
void mag_v_temp(double, double, double, int, int);


int main(){
  int rows, cols;
  double J_intra, J_inter, B, k1, k2, init_T;
  //lattice_t lattice;

  //rng = gsl_rng_alloc(gsl_rng_mt19937);
  //gsl_rng_set (rng, 0);

  clock_t begin = clock();

  initialize_lattice();
  M_v_B();
  //mag_v_temp(lattice, 5.0, .01, .05, 10000, 5000);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Execution time: %f \n", time_spent);

}

void initialize_lattice(){
  int i, j, k;

  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (rng, 0);
  //printf("In initialize_lattice\n");
  for(i = 0; i < ROWS; i++){
    for(j = 0; j < COLS; j++){
      gen_random_spin(&lattice[0][i][j]);
      gen_random_spin(&lattice[1][i][j]);
    }
  }
  //printf("end initialize_lattice\n");

}

void gen_random_spin(spin_t* spin){
    double x1, x2, mag_sq;
    gsl_rng_uniform(rng);
    //printf("In gen_random_spin\n");


    x1 = gsl_rng_uniform(rng)*2.0 - 1.0;
    x2 = gsl_rng_uniform(rng)*2.0 - 1.0;
    mag_sq = gsl_pow_2(x1) + gsl_pow_2(x2);
    //fprintf(stderr, "In gen_random_spin 2\n");


    while(mag_sq >= 1){
      x1 = gsl_rng_uniform(rng)*2.0 - 1.0;
      x2 = gsl_rng_uniform(rng)*2.0 - 1.0;
      mag_sq = gsl_pow_2(x1) + gsl_pow_2(x2);
    }
    //fprintf(stderr,"In gen_random_spin 3\n");

    spin->x = 2.0*x1*sqrt(1.0-mag_sq);
    spin->y = 2.0*x2*sqrt(1.0-mag_sq);
    spin->z = 1.0-2.0*mag_sq;

    //fprintf(stderr, "end gen_random_spin\n");

    //printf("x: %f\n", spin->x);
    //printf("y: %f\n", spin->y);
    //printf("z: %f\n", spin->z);

}

void simulate(int num_sweeps, double T){
  int i, num_accept = 0;
  for(i = 0; i < num_sweeps; i++){
    num_accept += sweep( T);
  }
  printf("Acceptance ratio: %f \n", num_accept/(num_sweeps*ROWS*COLS*2.0));
}

int sweep(double T){
  int i, num_accept = 0;
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


    /* Generate a perturbation of the spin at the chosen lattice site */
    perturb_spin(&temp_spin, &lattice[layer][row][col]);
    delta_E = calc_delta_E(&temp_spin, &lattice[layer][row][col], layer, row, col);
    random_num = gsl_rng_uniform(rng);
    //printf("random number: %f\n", random_num);
    if( !( (delta_E > 0) && (random_num >= gsl_sf_exp(-(1.0/T)*delta_E)) ) ){
      lattice[layer][row][col].x = temp_spin.x;
      lattice[layer][row][col].y = temp_spin.y;
      lattice[layer][row][col].z = temp_spin.z;
      num_accept += 1;
      //printf("Accepted. \n");
    }


  }

  return num_accept;

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

  //printf("Generated perturbation: %f, %f, %f \n", temp_spin->x, temp_spin->y, temp_spin->z);

  gsl_vector_free(spin_vector);
  gsl_vector_free(point_vector);
  gsl_vector_free(temp_vector);
  gsl_matrix_free(rot_matrix);
}

double calc_delta_E(spin_t* temp_spin, spin_t* spin, int layer, int row, int col){

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
  //printf("Delta E is %f \n", delta_E);

  gsl_vector_free(delta_vector);
  gsl_vector_free(neighbor_vector);
  gsl_vector_free(inter_vector);

  return delta_E;

}

void cool_lattice(double T){
    float curr_temp;
    int num_sweeps;
    curr_temp = INIT_T;
    while(curr_temp > T){
        curr_temp -= .02;
        printf("Cooling to %f \n", curr_temp);
        simulate(15000, curr_temp);
    }
}

double calc_magnetization(int layer){
      float mag, mag_spin;
      int i,j,k;
      mag = 0.0;
      mag_spin = 0.0;
      if(layer == -1){
        for(i=0; i < 2; i++)
            for(j=0; j < ROWS; j++)
                for(k = 0; k < COLS; k++)
                    mag += lattice[i][j][k].z;
        mag_spin = mag/(ROWS*COLS*2);
      }
      else{
        for(j=0; j < ROWS; j++)
            for(k = 0; k < COLS; k++)
                mag += lattice[layer][j][k].z;
        mag_spin = mag/(ROWS*COLS);
      }
      printf("Magnetization per spin (layer %d) is %f \n", layer, mag_spin);

      return mag_spin;
}

void M_v_B(){
    // Pass into M_v_B and mag_v_temp by reference
    double results[202][4];
    int sample_counter = 0;
    B_EXT = -1.0;
    cool_lattice(4.9);  // change to .1
    while(B_EXT < 1.0){
        printf("B: %f\n", B_EXT);
        simulate(15000, .1);
        // Measure magnetization
        results[sample_counter][0] = B_EXT;
        results[sample_counter][1] = calc_magnetization( -1);
        results[sample_counter][2] = calc_magnetization( 0);
        results[sample_counter][3] = calc_magnetization( 1);
        sample_counter += 1;

        B_EXT += .1;
    }

    while(B_EXT > -1.0){
      printf("B: %f\n", B_EXT);
        simulate(15000, .1);
        // Measure magnetization
        results[sample_counter][0] = B_EXT;
        results[sample_counter][1] = calc_magnetization( -1);
        results[sample_counter][2] = calc_magnetization( 0);
        results[sample_counter][3] = calc_magnetization( 1);
        sample_counter += 1;
        B_EXT -= .1;
    }
    
}

void mag_v_temp( double init_temp, double final_temp, double temp_step, int eq_time, int cor_time){

      double curr_temp = init_temp;
      double mag, mag1, mag2;

      while(curr_temp > final_temp){
          curr_temp -= temp_step;
          simulate( eq_time, curr_temp);

          simulate( cor_time, curr_temp);

          mag = calc_magnetization( -1);  // Mag
          mag1 = calc_magnetization( 0);  // Mag1
          mag2 = calc_magnetization( 1);  // Mag2

          printf("Current temp: %f \n", curr_temp);
          printf("Magnetization (total): %f \n", mag);
          printf("Magnetization (layer1): %f \n", mag1);
          printf("Magnetization (layer2): %f \n\n", mag2);
      }
      // **************Underflow error generated - check it out.
}
