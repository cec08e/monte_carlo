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

gcc -lgsl -lgslcblas heisenberg2d_1layer.c
ALT: gcc -fPIC -shared -o heisenberg2d_1layer.so -lgsl -lgslcblas heisenberg2d_1layer.c

************ N LAYER VERSION *****************
*/

#define ROWS 20       /* Number of rows in each lattice layer */
#define COLS 20       /* Number of columns in each lattice layer */
#define RADIUS .6     /* Radius of tangent disc in perturbing function */
#define INIT_T 5      /* Initial temperature */
#define DELTA_T .05   /* Annealing temp interval */
#define K 0           /* Anisotropy strength (negative for easy-axis) */
#define D .5          /* DM interaction strength */
#define NUM_L 1       /* Number of layers */

double B_EXT = -5.0;

int ANNEAL_TIME = 1000;                      /* Annealing speed */
int EQ_TIME = 1000;                          /* Number of equilibration sweeps */
int COR_TIME = 2000;                         /* Number of correlation sweeps */

typedef struct {
  double x;
  double y;
  double z;
} spin_t;             /* Spin type structure - x, y, and z components of spin */

typedef spin_t lattice_t[NUM_L][ROWS][COLS];

gsl_rng * rng;
lattice_t lattice;

float J_INTER[NUM_L];
float J_INTRA[NUM_L];
float K[NUM_L];

/* D vectors */
gsl_vector * D_n = gsl_vector_alloc(3); /* north neighbor */
gsl_vector * D_s = gsl_vector_alloc(3); /* south neighbor */
gsl_vector * D_e = gsl_vector_alloc(3); /* east neighbor */
gsl_vector * D_w = gsl_vector_alloc(3); /* west neighbor */

void initialize_lattice();
void initialize_params();
void gen_random_spin(spin_t*);
void simulate( int, double);
void perturb_spin(spin_t* , spin_t*);
double calc_delta_E(spin_t* , spin_t*, int, int, int);
void cross_product(const gsl_vector*, const gsl_vector*, gsl_vector*);
void cool_lattice(double);
double calc_magnetization(int);
int M_v_B(double**);

void initialize_lattice(){
  int i, j, k;

  initialize_params();

  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (rng, time(NULL));
  for(l = 0; l < NUM_L; l++){
    for(i = 0; i < ROWS; i++){
      for(j = 0; j < COLS; j++){
        gen_random_spin(&lattice[l][i][j]);
      }
    }
  }
}

void initialize_params(){
  /* Initialize D vector values */

  init_D_vec(D_n, 0, D, 0);
  init_D_vec(D_s, 0, -D, 0);
  init_D_vec(D_e, D, 0, 0);
  init_D_vec(D_w, -D, 0, 0);


  /* Change K, J_inter and J_intra parameters here */
  int j;
  for(j = 0; j < NUM_L; j++){
    K[j] = .05;
    J_INTER[j] = .1;
    J_INTRA[j] = 1.0;
  }

  J_INTER[NUM_L-1] = 0; /* No interaction between 1st and last layer */

  /* Example, introducing small increased anisotropy on top layer:

  K[0] = .08

  */


}

void init_D_vec(gsl_vector* D_vec, x, y, z){

  gsl_vector_set(D_vec, 0, x);
  gsl_vector_set(D_vec, 1, y);
  gsl_vector_set(D_vec, 2, z);


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

}

void simulate(int num_sweeps, double T){
  int i, num_accept = 0;
  for(i = 0; i < num_sweeps; i++){
    num_accept += sweep(T);
  }
  /*
  if(num_sweeps > 1){
    printf("Acceptance ratio: %f \n", num_accept/(num_sweeps*ROWS*COLS*NUM_L));
  }
  */
}

int sweep(double T){
  int i, num_accept = 0;
  int layer, row, col;
  double delta_E;
  double random_num;
  spin_t temp_spin;

  /* Perform as many MC steps as there are lattice sites */
  for(i=0; i < NUM_L*ROWS*COLS; i++){

    /* Choose a random spin site on a random lattice */
    layer = gsl_rng_get(rng) % NUM_L;   /******************************************************/
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
  double delta_D;
  int i;

  /* FIRST TERM */
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

  gsl_blas_ddot(delta_vector, neighbor_vector, &delta_dot_neighbor);
  /* END FIRST TERM */

  /* SECOND TERM */

  gsl_vector* inter_vector = gsl_vector_alloc(3);

  gsl_vector_set(inter_vector, 0, J_INTER[layer]*lattice[(layer+1)%NUM_L][row][col].x + J_INTER[(((layer-1)%NUM_L)+NUM_L)%NUM_L]*lattice[(((layer-1)%NUM_L)+NUM_L)%NUM_L][row][col].x);
  gsl_vector_set(inter_vector, 1, J_INTER[layer]*lattice[(layer+1)%NUM_L][row][col].y + J_INTER[(((layer-1)%NUM_L)+NUM_L)%NUM_L]*lattice[(((layer-1)%NUM_L)+NUM_L)%NUM_L][row][col].y);
  gsl_vector_set(inter_vector, 2, J_INTER[layer]*lattice[(layer+1)%NUM_L][row][col].z + J_INTER[(((layer-1)%NUM_L)+NUM_L)%NUM_L]*lattice[(((layer-1)%NUM_L)+NUM_L)%NUM_L ][row][col].z);

  gsl_blas_ddot(delta_vector, inter_vector, &delta_dot_inter);
  /* END SECOND TERM */

  /* THIRD TERM */
  /* Calculate anisotropy change */
  delta_a = K[layer]*(gsl_pow_2(temp_spin->z) - gsl_pow_2(spin->z));
  /* END THIRD TERM */

  /* FOURTH TERM */
  /* TO DO: DELTA D CALCULATION */
  delta_D = 0;
  /*
  for (i = 0; i < 3; i++){
    printf ("neighbor_vector_%d = %g\n", i, gsl_vector_get (neighbor_vector, i));
  }
  printf("delta_dot_neighbor: %f\n", delta_dot_neighbor);
  printf("delta_dot_inter: %f\n", delta_dot_inter);
  printf("delta_a: %f\n", delta_a);
  */

  /* delta_D = dot((0,self.D,0), cross(delta_spin, self.lattice[(row-1)%self.rows][col]))   # north neighbor
      delta_D += dot((0,-1*self.D,0), cross(delta_spin, self.lattice[(row+1)%self.rows][col]))    # south neighbor
      delta_D += dot((-1*self.D,0,0), cross(delta_spin, self.lattice[row][(col-1)%self.cols]))    # west neighbor
      delta_D += dot((self.D,0,0), cross(delta_spin, self.lattice[row][(col+1)%self.cols]))    # east neighbor */
  //gsl_blas_ddot(delta_vector, neighbor_vector, &delta_D);

  /* END FOURTH TERM */




  //delta_E = -J_INTRA*delta_dot_neighbor + J_INTER*delta_dot_inter + delta_a - B_EXT*gsl_vector_get(delta_vector,2);
  delta_E = -J_INTRA[layer]*delta_dot_neighbor + delta_dot_inter + delta_a + delta_D - B_EXT*gsl_vector_get(delta_vector,2);

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
      curr_temp -= DELTA_T;
      simulate(ANNEAL_TIME, curr_temp);
  }

}

double calc_magnetization(int layer){
      float mag, mag_spin;
      int i,j,k;
      mag = 0.0;
      mag_spin = 0.0;
      if(layer == -1){
        for(i=0; i < NUM_L; i++)
            for(j=0; j < ROWS; j++)
                for(k = 0; k < COLS; k++)
                    mag += lattice[i][j][k].z;
        mag_spin = mag/(ROWS*COLS*NUM_L);
      }
      else{
        for(j=0; j < ROWS; j++)
            for(k = 0; k < COLS; k++)
                mag += lattice[layer][j][k].z;
        mag_spin = mag/(ROWS*COLS);
      }
      return mag_spin;
}


/* EXPERIMENTS */
int M_v_B(double** results){
    int cor_count = 0;
    //clock_t begin = clock();
    //FILE *f = fopen(SIM_CONFIG, "a");
    //fprintf(f, "Simulation %d: Size = %d, J_inter = {%f,%f,%f,%f}, K = {%f,%f,%f,%f}, J_intra = {%f,%f,%f,%f}, T=.15, Steps=5000, dB = .004\n", SIM_NUM, ROWS, J_INTER[0], J_INTER[1],
    //J_INTER[2], J_INTER[3], K[0], K[1], K[2], K[3], J_INTRA[0], J_INTRA[1], J_INTRA[2], J_INTRA[3]);
    //fclose(f);
    int sample_counter = 0;
    int i;

    B_EXT = -.2;
    delta_B = .005;

    cool_lattice(.15);
    while(B_EXT < .2){
        printf("B: %f\n", B_EXT);
        simulate(EQ_TIME, .15);
        // Measure magnetization
        results[sample_counter][0] = B_EXT;
        for(i=0; i <= NUM_L; i++)
          results[sample_counter][i+1] = calc_magnetization(i-1);

        // Take average over 1000 sweeps
        for(cor_count = 1; cor_count < COR_TIME; cor_count++){
          simulate(1, .15);
          for(i=0; i <= NUM_L; i++)
            results[sample_counter][i+1] = calc_magnetization(i-1);
        }
        for(i=0; i <= NUM_L; i++)
          results[sample_counter][i+1] += results[sample_counter][i+1]/COR_TIME;

        sample_counter += 1;

        B_EXT += delta_B;
    }

    while(B_EXT > -.2){
        printf("B: %f\n", B_EXT);
        simulate(EQ_TIME, .15);
        // Measure magnetization
        results[sample_counter][0] = B_EXT;
        for(i=0; i <= NUM_L; i++)
          results[sample_counter][i+1] = calc_magnetization(i-1);

        // Take average over 1000 sweeps
        for(cor_count = 1; cor_count < COR_TIME; cor_count++){
          simulate(1, .15);
          for(i=0; i <= NUM_L; i++)
            results[sample_counter][i+1] = calc_magnetization(i-1);
        }
        for(i=0; i <= NUM_L; i++)
          results[sample_counter][i+1] += results[sample_counter][i+1]/COR_TIME;
        sample_counter += 1;
        B_EXT -= delta_B;
    }

    //clock_t end = clock();
    //double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("Execution time: %f \n", time_spent);

    return sample_counter;

}

/* Helper function for energy calculation */

void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product)
{
        double p1 = gsl_vector_get(u, 1)*gsl_vector_get(v, 2)
                - gsl_vector_get(u, 2)*gsl_vector_get(v, 1);

        double p2 = gsl_vector_get(u, 2)*gsl_vector_get(v, 0)
                - gsl_vector_get(u, 0)*gsl_vector_get(v, 2);

        double p3 = gsl_vector_get(u, 0)*gsl_vector_get(v, 1)
                - gsl_vector_get(u, 1)*gsl_vector_get(v, 0);

        gsl_vector_set(product, 0, p1);
        gsl_vector_set(product, 1, p2);
        gsl_vector_set(product, 2, p3);
}
