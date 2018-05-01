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

************ 5 LAYER VERSION *****************
*/

#define ROWS 10       /* Number of rows in each lattice layer */
#define COLS 10       /* Number of columns in each lattice layer */
#define RADIUS .6   /* Radius of tangent disc in perturbing function */
#define J_INTRA 1
/*
#define J_INTER_12 .1
#define J_INTER_23 .1
#define J_INTER_34 .1
*/
double J_INTER[4] = {.1, .1, .1, .1, 0};

// Same K on all layers - try .2
// Down to small interlayer
// Same configuration sweeping up and back
// 5 layers


#define INIT_T 5

// Play around with inter, k , etc...- keep intra at 1
// Plot neel mag vs temp to find neel temp
//


double B_EXT = -5;
/*
double K1 = -1;
double K2 = -1;
double K3 = -1;
double K4 = -1;
*/
double K[4] = {-.2, -.2, -.2, -.2, -.2};
int EQ_TIME = 5000;
int COR_TIME = 5000;

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

typedef spin_t lattice_t[5][ROWS][COLS];

gsl_rng * rng;
lattice_t lattice;


void initialize_lattice();
void gen_random_spin(spin_t*);
void simulate( int, double);
int sweep( double);
void perturb_spin(spin_t* , spin_t*);
double calc_delta_E( spin_t* , spin_t*, int, int, int);
double calc_magnetization( int);
int M_v_B(double**);
int M_v_K(double**);
int M_v_T(double**, double, double, double);


int main(){
  int rows, cols;
  //double J_intra, J_inter, B, k1, k2, init_T;
  //lattice_t lattice;

  //rng = gsl_rng_alloc(gsl_rng_mt19937);
  //gsl_rng_set (rng, 0);

  clock_t begin = clock();

  initialize_lattice();
  //M_v_B();
  //mag_v_temp(lattice, 5.0, .01, .05, 10000, 5000);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  //printf("Execution time: %f \n", time_spent);

}

void initialize_lattice(){
  int i, j, k;

  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (rng, 0);
  for(i = 0; i < ROWS; i++){
    for(j = 0; j < COLS; j++){
      gen_random_spin(&lattice[0][i][j]);
      gen_random_spin(&lattice[1][i][j]);
      gen_random_spin(&lattice[2][i][j]);
      gen_random_spin(&lattice[3][i][j]);
      gen_random_spin(&lattice[4][i][j]);

    }
  }
}

void gen_random_spin(spin_t* spin){
    double x1, x2, mag_sq;
    gsl_rng_uniform(rng);

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

    /*
    printf("x: %f\n", spin->x);
    printf("y: %f\n", spin->y);
    printf("z: %f\n", spin->z);
    printf("Magnitude: %f\n\n", sqrt(gsl_pow_2(spin->x)
                + gsl_pow_2(spin->y)
                + gsl_pow_2(spin->z)));
    */

}

void simulate(int num_sweeps, double T){
  int i, num_accept = 0;
  for(i = 0; i < num_sweeps; i++){
    num_accept += sweep(T);
  }
  printf("Acceptance ratio: %f \n", num_accept/(num_sweeps*ROWS*COLS*5.0));
}

int sweep(double T){
  int i, num_accept = 0;
  int layer, row, col;
  double delta_E;
  double random_num;
  spin_t temp_spin;


  /* Perform as many MC steps as there are lattice sites */
  for(i=0; i < 5*ROWS*COLS; i++){

    /* Choose a random spin site on a random lattice */
    layer = gsl_rng_uniform_int(rng,5);
    row = gsl_rng_get(rng) % ROWS;
    col = gsl_rng_get(rng) % COLS;

    //printf("Chosen site: %d, %d, %d \n", layer, row, col);
    //printf("Spin is: (%f, %f, %f)\n", lattice[layer][row][col].x, lattice[layer][row][col].y, lattice[layer][row][col].z);


    /* Generate a perturbation of the spin at the chosen lattice site */
    perturb_spin(&temp_spin, &lattice[layer][row][col]);
    delta_E = calc_delta_E(&temp_spin, &lattice[layer][row][col], layer, row, col);
    random_num = gsl_rng_uniform(rng);
    //printf("Delta E is %f\n", delta_E);
    //printf("Random number: %f\n", random_num);
    //printf("Exponential: %f \n", gsl_sf_exp(-(1.0/T)*delta_E));
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
  //printf("Delta vector is (%f,%f,%f) \n", gsl_vector_get(delta_vector,0), gsl_vector_get(delta_vector,1), gsl_vector_get(delta_vector,2));
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
  //printf("Jinter dot with layer %d \n", (layer+1)%2);

  //// ADD INTER LAYER calc changes

  gsl_vector_set(inter_vector, 0, J_INTER[layer]*lattice[(layer+1)%4][row][col].x + J_INTER[(((layer-1)%4)+4)%4]*lattice[(((layer-1)%4)+4)%4][row][col].x);
  gsl_vector_set(inter_vector, 1, J_INTER[layer]*lattice[(layer+1)%4][row][col].y + J_INTER[(((layer-1)%4)+4)%4]*lattice[(((layer-1)%4)+4)%4][row][col].y);
  gsl_vector_set(inter_vector, 2, J_INTER[layer]*lattice[(layer+1)%4][row][col].z + J_INTER[(((layer-1)%4)+4)%4]*lattice[(((layer-1)%4)+4)%4][row][col].z);

  gsl_blas_ddot(delta_vector, inter_vector, &delta_dot_inter);


  /* Calculate anisotropy change */
  delta_a = K[layer]*(gsl_pow_2(temp_spin->z) - gsl_pow_2(spin->z));

  /*
  for (i = 0; i < 3; i++){
    printf ("neighbor_vector_%d = %g\n", i, gsl_vector_get (neighbor_vector, i));
  }
  printf("delta_dot_neighbor: %f\n", delta_dot_neighbor);
  printf("delta_dot_inter: %f\n", delta_dot_inter);
  printf("delta_a: %f\n", delta_a);
  */

  //delta_E = -J_INTRA*delta_dot_neighbor + J_INTER*delta_dot_inter + delta_a - B_EXT*gsl_vector_get(delta_vector,2);
  delta_E = -J_INTRA*delta_dot_neighbor + delta_dot_inter + delta_a - B_EXT*gsl_vector_get(delta_vector,2);

  //printf("Delta E is %f \n", delta_E);

  gsl_vector_free(delta_vector);
  gsl_vector_free(neighbor_vector);
  gsl_vector_free(inter_vector);

  return delta_E;

}

void cool_lattice(double T){
    float curr_temp;
    curr_temp = INIT_T;
    while(curr_temp > T){
        curr_temp -= .035;
        printf("Cooling to %f \n", curr_temp);
        simulate(2000, curr_temp);
    }
}

double calc_magnetization(int layer){
      float mag, mag_spin;
      int i,j,k;
      mag = 0.0;
      mag_spin = 0.0;
      if(layer == -1){
        for(i=0; i < 5; i++)
            for(j=0; j < ROWS; j++)
                for(k = 0; k < COLS; k++)
                    mag += lattice[i][j][k].z;
        mag_spin = mag/(ROWS*COLS*4);
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

int M_v_B(double** results){

    int sample_counter = 0;
    B_EXT = -.7;
    cool_lattice(.15);
    while(B_EXT < .7){
        printf("B: %f\n", B_EXT);
        simulate(10000, .15);
        // Measure magnetization
        results[sample_counter][0] = B_EXT;
        results[sample_counter][1] = calc_magnetization( -1);
        results[sample_counter][2] = calc_magnetization( 0);
        results[sample_counter][3] = calc_magnetization( 1);
        results[sample_counter][4] = calc_magnetization( 2);
        results[sample_counter][5] = calc_magnetization( 3);
        results[sample_counter][6] = calc_magnetization( 4);

        sample_counter += 1;

        B_EXT += .003;
    }

    while(B_EXT > -.7){
        printf("B: %f\n", B_EXT);
        simulate(10000, .15);
        // Measure magnetization
        results[sample_counter][0] = B_EXT;
        results[sample_counter][1] = calc_magnetization( -1);
        results[sample_counter][2] = calc_magnetization( 0);
        results[sample_counter][3] = calc_magnetization( 1);
        results[sample_counter][4] = calc_magnetization( 2);
        results[sample_counter][5] = calc_magnetization( 3);
        results[sample_counter][6] = calc_magnetization( 4);

        sample_counter += 1;
        B_EXT -= .003;
    }

    return sample_counter;

}

int M_v_K(double** results){

    int sample_counter = 0;
    //B_EXT = .6;   /* Critical Switching Field */
    double K_val = -1.0;
    cool_lattice(.1);
    while(K_val < 0){
        K[0] = K_val;
        K[1] = K_val;
        K[2] = K_val;
        K[3] = K_val;
        printf("K: %f\n", K_val);
        simulate(10000, .1);
        // Measure magnetization
        results[sample_counter][0] = K_val;
        results[sample_counter][1] = calc_magnetization( -1);
        results[sample_counter][2] = calc_magnetization( 0);
        results[sample_counter][3] = calc_magnetization( 1);
        results[sample_counter][4] = calc_magnetization( 2);
        results[sample_counter][5] = calc_magnetization( 3);
        results[sample_counter][6] = calc_magnetization( 4);

        sample_counter += 1;

        K_val += .005;
    }

    while(K_val > -1.0){
        K[0] = K_val;
        K[1] = K_val;
        K[2] = K_val;
        K[3] = K_val;
        printf("K: %f\n", K_val);
        simulate(10000, .1);
        // Measure magnetization
        results[sample_counter][0] = K_val;
        results[sample_counter][1] = calc_magnetization( -1);
        results[sample_counter][2] = calc_magnetization( 0);
        results[sample_counter][3] = calc_magnetization( 1);
        results[sample_counter][4] = calc_magnetization( 2);
        results[sample_counter][5] = calc_magnetization( 3);
        results[sample_counter][6] = calc_magnetization( 4);

        sample_counter += 1;
        K_val -= .005;
    }

    return sample_counter;

}

/*******************************************************************************/
//   The M_v_T function generates a series of magnetization
//   measurements for a series of time steps.
// UPDATE for 5 Layer
/*******************************************************************************/

int M_v_T(double** results, double init_temp, double final_temp, double temp_step){
      double curr_temp = init_temp;
      double mag, mag1, mag2, neel_mag;
      int sample_counter = 0;
      printf("init_temp: %f\n", init_temp);
      printf("final_temp: %f\n", final_temp);
      printf("temp_step: %f\n", temp_step);

      while(curr_temp > final_temp){
          printf("Current temp: %f\n", curr_temp);

          simulate(EQ_TIME, curr_temp);
          simulate(COR_TIME, curr_temp);

          results[sample_counter][0] = curr_temp;
          results[sample_counter][1] = calc_magnetization(-1);  // Mag
          results[sample_counter][2] = calc_magnetization(0);  // Mag1
          results[sample_counter][3] = calc_magnetization(1);  // Mag2
          results[sample_counter][4] = results[sample_counter][2] - results[sample_counter][3] ;
          /*
          printf("Current temp: %f \n", curr_temp);
          printf("Magnetization (total): %f \n", mag);
          printf("Magnetization (layer1): %f \n", mag1);
          printf("Magnetization (layer2): %f \n\n", mag2);
          */
          curr_temp -= temp_step;
          sample_counter += 1;
      }
      // **************Underflow error generated - check it out.
      return sample_counter;
}