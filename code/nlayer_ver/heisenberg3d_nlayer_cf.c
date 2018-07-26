#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <time.h>

/*
Requirements: GNU Scientific Library (gsl) and CBLAS

gcc -lgsl -lgslcblas heisenberg2d_1layer.c
ALT: gcc -fPIC -shared -o heisenberg2d_1layer.so -lgsl -lgslcblas heisenberg2d_1layer.c

************ N LAYER VERSION *****************
*/

/*********************************************/
/*    Simulation Parameters                  */
/*    May be set here or passed through      */
/*    initialize lattice.                    */
/*********************************************/

int SIM_NUM = 185;    /* Simulation number - appears in both configuration and result files */

int NUM_L = 1;           /* Number of layers */
int ROWS = 20;        /* Number of rows in each lattice layer */
int COLS = 20;        /* Number of columns in each lattice layer */

double RADIUS = .6;       /* Radius of tangent disc in perturbing function */
double INIT_T = 4.0;      /* Initial temperature */
double TEMP = .5;         /* Final temperature */
double DELTA_T = .1;      /* Annealing temp interval */
double DELTA_B = .004;    /* B sweeping speed */
double D = 0.3;           /* DM interaction strength */
double B_EXT = .2;

float J_INTER[NUM_L];
float J_INTRA[NUM_L];
float K[NUM_L];

int OVER_FLAG = 2;        /* Number of OR sweeps performed after every MC sweep */

int ANNEAL_TIME = 2000;         /* Temperature annealing time */
int EQ_TIME = 100000;           /* Number of equilibration sweeps */
int COR_TIME = 100;             /* Number of correlation sweeps */

/********************************************/
/*    End simulation parameters.            */
/********************************************/

#define SIM_CONFIG "sim_results/sim_config.txt"

typedef struct {
  double x;
  double y;
  double z;
} spin_t;             /* Spin type structure - x, y, and z components of spin */

typedef spin_t lattice_t[NUM_L][ROWS][COLS];

gsl_rng * rng;
lattice_t lattice;
lattice_t lattice_copy;

/* D vectors */
gsl_vector * D_vec[4];

void initialize_lattice(char *);
void parse_config_file(char *);
void initialize_params();
void init_D_vec();
void gen_random_spin(spin_t*);
void simulate( int, double);
int sweep(double T);
void perturb_spin(spin_t* , spin_t*);
double calc_delta_E(spin_t* , spin_t*, int, int, int);
void cross_product(const gsl_vector*, const gsl_vector*, gsl_vector*);
void cool_lattice(double);
void record_lattice(spin_t***);
void overrelax();
void eff_project(spin_t*, int , int, int);
double calc_magnetization(int);
int M_v_B(double**);
double calc_TC();
double calc_solid_angle(spin_t n1, spin_t n2, spin_t n3);
void sample_mag(double*, int);
int mag_v_temp(double**);
int suscept_v_temp(double**);
void test_lattice_TC();


int main(){
  printf("Testing");
  test_lattice_TC();
}

void initialize_lattice(char * config_fn){
  int l, i, j;
  for(i=0; i < 4; i++)
    D_vec[i] = gsl_vector_alloc(3

  init_D_vec();

  if(!strcmp(config_fn, "")){
    printf("Using configuration file: %s\n", config_fn);
    parse_config_file(config_fn);
  }

  //initialize_params();

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

void parse_config_file(char * config_fn){
  /* Format of config file:

  */
}

void initialize_params(){
  /* Initialize D vector values */

  init_D_vec();
  //init_D_vec(D_n, 0, D, 0);
  //init_D_vec(D_s, 0, -D, 0);
  //init_D_vec(D_e, D, 0, 0);
  //init_D_vec(D_w, -D, 0, 0);


  /* Change K, J_inter and J_intra parameters here */
  int j;
  for(j = 0; j < NUM_L; j++){
    K[j] = -0;
    J_INTER[j] = .1;
    J_INTRA[j] = 1.0;
  }

  //J_INTER[0] = .1;   /* TOP LAYERS */
  J_INTER[NUM_L-1] = 0; /* No interaction between 1st and last layer */

  /* Example, introducing small increased anisotropy on top layer:

  K[0] = .08

  */


}

void init_D_vec(){

  /*
  gsl_vector_set(D_vec[0], 0, 0);
  gsl_vector_set(D_vec[0], 1, D);
  gsl_vector_set(D_vec[0], 2, 0);


  gsl_vector_set(D_vec[1], 0, 0);
  gsl_vector_set(D_vec[1], 1, -D);
  gsl_vector_set(D_vec[1], 2, 0);


  gsl_vector_set(D_vec[2], 0, -D);
  gsl_vector_set(D_vec[2], 1, 0);
  gsl_vector_set(D_vec[2], 2, 0);


  gsl_vector_set(D_vec[3], 0, D);
  gsl_vector_set(D_vec[3], 1, 0);
  gsl_vector_set(D_vec[3], 2, 0);
  */

  gsl_vector_set(D_vec[0], 0, D);
  gsl_vector_set(D_vec[0], 1, 0);
  gsl_vector_set(D_vec[0], 2, 0);


  gsl_vector_set(D_vec[1], 0, -D);
  gsl_vector_set(D_vec[1], 1, 0);
  gsl_vector_set(D_vec[1], 2, 0);


  gsl_vector_set(D_vec[2], 0, 0);
  gsl_vector_set(D_vec[2], 1, -D);
  gsl_vector_set(D_vec[2], 2, 0);


  gsl_vector_set(D_vec[3], 0, 0);
  gsl_vector_set(D_vec[3], 1, D);
  gsl_vector_set(D_vec[3], 2, 0);


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
  int i, num_or, num_accept = 0;
  for(i = 0; i < num_sweeps; i++){
    //printf("Sweep #: %d\n", i);
    num_accept += sweep(T);
    for(num_or = 0; num_or < OVER_FLAG; num_or++){
      //printf("Overrelaxing.");
      overrelax();
    }
  }

  if(num_sweeps > 1){
    printf("Acceptance ratio: %f \n", num_accept/(num_sweeps*ROWS*COLS*NUM_L*1.0));
  }

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

  gsl_vector* neighbors[4];
  gsl_vector* cross_temp = gsl_vector_alloc(3);
  double delta_dot_neighbor;
  double delta_dot_inter;
  double delta_a;
  double delta_E;
  double delta_D;
  double temp;
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

  /* NORTH */
  neighbors[0] = gsl_vector_alloc(3);
  gsl_vector_set(neighbors[0], 0, lattice[layer][(((row-1)%ROWS) + ROWS) % ROWS][col].x);
  gsl_vector_set(neighbors[0], 1, lattice[layer][(((row-1)%ROWS) + ROWS) % ROWS][col].y);
  gsl_vector_set(neighbors[0], 2, lattice[layer][(((row-1)%ROWS) + ROWS) % ROWS][col].z);

  /* SOUTH */
  neighbors[1] = gsl_vector_alloc(3);
  gsl_vector_set(neighbors[1], 0, lattice[layer][(((row+1)%ROWS) + ROWS) % ROWS][col].x);
  gsl_vector_set(neighbors[1], 1, lattice[layer][(((row+1)%ROWS) + ROWS) % ROWS][col].y);
  gsl_vector_set(neighbors[1], 2, lattice[layer][(((row+1)%ROWS) + ROWS) % ROWS][col].z);

  /* EAST */
  neighbors[2] = gsl_vector_alloc(3);
  gsl_vector_set(neighbors[2], 0, lattice[layer][row][(((col+1)%COLS) + COLS) % COLS].x);
  gsl_vector_set(neighbors[2], 1, lattice[layer][row][(((col+1)%COLS) + COLS) % COLS].y);
  gsl_vector_set(neighbors[2], 2, lattice[layer][row][(((col+1)%COLS) + COLS) % COLS].z);

  neighbors[3] = gsl_vector_alloc(3);
  gsl_vector_set(neighbors[3], 0, lattice[layer][row][(((col-1)%COLS) + COLS) % COLS].x);
  gsl_vector_set(neighbors[3], 1, lattice[layer][row][(((col-1)%COLS) + COLS) % COLS].y);
  gsl_vector_set(neighbors[3], 2, lattice[layer][row][(((col-1)%COLS) + COLS) % COLS].z);

  delta_D = 0;

  for(i = 0; i < 4; i++){
    cross_product(delta_vector, neighbors[i], cross_temp);
    gsl_blas_ddot(D_vec[i], cross_temp, &temp);
    delta_D += temp;
  }

  //printf("Delta D is %f \n", delta_D);
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


  //printf("Inter term: %f \n", delta_dot_inter);
  //printf("Anisotropy term: %f \n", delta_a);
  //printf("DM term: -%f \n", delta_D);

  //delta_E = -J_INTRA*delta_dot_neighbor + J_INTER*delta_dot_inter + delta_a - B_EXT*gsl_vector_get(delta_vector,2);
  delta_E = -J_INTRA[layer]*delta_dot_neighbor + delta_dot_inter + delta_a + delta_D - B_EXT*gsl_vector_get(delta_vector,2);

  //printf("Delta E is %f \n", delta_E);

  gsl_vector_free(delta_vector);
  gsl_vector_free(neighbor_vector);
  gsl_vector_free(inter_vector);
  gsl_vector_free(cross_temp);
  gsl_vector_free(neighbors[0]);
  gsl_vector_free(neighbors[1]);
  gsl_vector_free(neighbors[2]);
  gsl_vector_free(neighbors[3]);

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

void overrelax(){
  /* Overrelaxation step performed between every sweep.

    Reflect every spin on the lattice to the other side of the effective field vector.

    QUESTION: Does overrelaxing need to be performed using the entire old-state? Jiadong
    uses already overrelaxed spin configurations in future overrelaxation procedures.

  */
  spin_t temp;
  int i, j, k;
  for(i = 0; i < NUM_L; i++)
    for(j=0; j < ROWS; j++)
      for(k = 0; k < COLS; k++){
        eff_project(&temp, i, j, k);
        /* Reflect spin (j,k) on layer i */
        lattice[i][j][k].x = 2*temp.x - lattice[i][j][k].x;
        lattice[i][j][k].y = 2*temp.y - lattice[i][j][k].y;
        lattice[i][j][k].z = 2*temp.z - lattice[i][j][k].z;
      }
  /*
  for(i = 0; i < NUM_L; i++)
    for(j=0; j < ROWS; j++)
      for(k = 0; k < COLS; k++){
            lattice[i][j][k].x = lattice_copy[i][j][k].x;
            lattice[i][j][k].y = lattice_copy[i][j][k].y;
            lattice[i][j][k].z = lattice_copy[i][j][k].z;
      }

  */

}

void eff_project(spin_t* field, int layer, int row, int col){

  /* Calculate and store V(s*V/|V|^2) for spin at (j,k) on layer i */

  gsl_vector * V = gsl_vector_alloc(3);
  gsl_vector * s = gsl_vector_alloc(3);
  //gsl_vector * N = gsl_vector_alloc(3);
  //gsl_vector * S = gsl_vector_alloc(3);
  //gsl_vector * E = gsl_vector_alloc(3);
  //gsl_vector * W = gsl_vector_alloc(3);
  gsl_vector * cross_temp = gsl_vector_alloc(3);
  int i;
  double V_x = 0, V_y = 0, V_z = 0;
  double temp;

  gsl_vector* neighbors[4];

  gsl_vector_set(s, 0, lattice[layer][row][col].x);
  gsl_vector_set(s, 1, lattice[layer][row][col].y);
  gsl_vector_set(s, 2, lattice[layer][row][col].z);


  /* NORTH */
  neighbors[0] = gsl_vector_alloc(3);
  gsl_vector_set(neighbors[0], 0, lattice[layer][(((row-1)%ROWS) + ROWS) % ROWS][col].x);
  gsl_vector_set(neighbors[0], 1, lattice[layer][(((row-1)%ROWS) + ROWS) % ROWS][col].y);
  gsl_vector_set(neighbors[0], 2, lattice[layer][(((row-1)%ROWS) + ROWS) % ROWS][col].z);

  /* SOUTH */
  neighbors[1] = gsl_vector_alloc(3);
  gsl_vector_set(neighbors[1], 0, lattice[layer][(((row+1)%ROWS) + ROWS) % ROWS][col].x);
  gsl_vector_set(neighbors[1], 1, lattice[layer][(((row+1)%ROWS) + ROWS) % ROWS][col].y);
  gsl_vector_set(neighbors[1], 2, lattice[layer][(((row+1)%ROWS) + ROWS) % ROWS][col].z);

  /* EAST */
  neighbors[2] = gsl_vector_alloc(3);
  gsl_vector_set(neighbors[2], 0, lattice[layer][row][(((col+1)%COLS) + COLS) % COLS].x);
  gsl_vector_set(neighbors[2], 1, lattice[layer][row][(((col+1)%COLS) + COLS) % COLS].y);
  gsl_vector_set(neighbors[2], 2, lattice[layer][row][(((col+1)%COLS) + COLS) % COLS].z);

  /* WEST */
  neighbors[3] = gsl_vector_alloc(3);
  gsl_vector_set(neighbors[3], 0, lattice[layer][row][(((col-1)%COLS) + COLS) % COLS].x);
  gsl_vector_set(neighbors[3], 1, lattice[layer][row][(((col-1)%COLS) + COLS) % COLS].y);
  gsl_vector_set(neighbors[3], 2, lattice[layer][row][(((col-1)%COLS) + COLS) % COLS].z);



  /* Add together neighbors with J_intra, and cross products of neighbors with D, and Bz*/

  for(i = 0; i < 4; i++){
    V_x += -J_INTRA[layer]*gsl_vector_get(neighbors[i], 0);
    V_y += -J_INTRA[layer]*gsl_vector_get(neighbors[i], 1);
    V_z += -J_INTRA[layer]*gsl_vector_get(neighbors[i], 2);

    cross_product(neighbors[i], D_vec[i], cross_temp);
    V_x += gsl_vector_get(cross_temp, 0);
    V_y += gsl_vector_get(cross_temp, 1);
    V_z += gsl_vector_get(cross_temp, 2);

  }
  V_z += -B_EXT;


  gsl_vector_set(V, 0, V_x);
  gsl_vector_set(V, 1, V_y);
  gsl_vector_set(V, 2, V_z);

  gsl_blas_ddot(s, V, &temp);
  temp = temp/(gsl_pow_2(V_x) + gsl_pow_2(V_y) + gsl_pow_2(V_z));

  field->x = V_x*temp;
  field->y = V_y*temp;
  field->z = V_z*temp;

  gsl_vector_free(V);
  gsl_vector_free(s);
  gsl_vector_free(cross_temp);
  gsl_vector_free(neighbors[0]);
  gsl_vector_free(neighbors[1]);
  gsl_vector_free(neighbors[2]);
  gsl_vector_free(neighbors[3]);
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

void record_lattice(spin_t*** record){
  int i, j, k, n;
  double TC = 0;


  /* Cool the lattice to T  and record the lattice configuration after 5000 sweeps */

  cool_lattice(TEMP);
  simulate(EQ_TIME, TEMP);
  for(i = 0; i < 3000; i++){
    simulate(COR_TIME, TEMP);
    TC += calc_TC();
    printf("TC average %d: %f\n", i+1, TC/(i+1));
  }
  TC = TC/3000.0;

  /* Record lattice configuration to results */
  for(i = 0; i < NUM_L; i++)
    for(j = 0; j < ROWS; j++)
      for(k = 0; k < COLS; k++){
        record[i][j][k].x = lattice[i][j][k].x;
        record[i][j][k].y = lattice[i][j][k].y;
        record[i][j][k].z = lattice[i][j][k].z;
      }

  FILE *f = fopen(SIM_CONFIG, "a");
  fprintf(f, "Simulation %d\n", SIM_NUM);
  fprintf(f, "Lattice Record %d: Size = %d\n", SIM_NUM, ROWS);
  fprintf(f,"\tD = %f\n", D);
  fprintf(f,"\tB = %f\n", B_EXT);
  fprintf(f,"\tT = %f\n", TEMP);
  fprintf(f,"\tOverrelaxation: %d\n", OVER_FLAG);
  fprintf(f,"\tCOR time = %d sweeps \n", COR_TIME);
  fprintf(f,"\tEQ time = %d sweeps \n", EQ_TIME);
  fprintf(f,"\tTC averaged over 100 correlation times.\n");

  for(n = 0; n < NUM_L; n++)
    fprintf(f, "\tJ_AF[%d] = %f\n\tK[%d] = %f\n", n, J_INTER[n], n, K[n]);
  fprintf(f, "\tTopological charge: %f \n", TC);
  fclose(f);

  printf("Topological charge is %f \n", TC);

}


/* EXPERIMENTS */
int M_v_B(double** results){
    int cor_count = 0;
    int n = 0;
    //clock_t begin = clock();
    FILE *f = fopen(SIM_CONFIG, "a");
    fprintf(f, "Simulation %d: Size = %d, T=.15, Steps=%d, dB = %f\n", SIM_NUM, ROWS, COR_TIME, DELTA_B);
    for(n = 0; n < NUM_L; n++)
      fprintf(f, "\tJ_AF[%d] = %f\n\tK[%d] = %f\n", n, J_INTER[n], n, K[n]);
    fclose(f);
    int sample_counter = 0;
    int i;

    B_EXT = -.2;
    //float delta_B = .004;

    cool_lattice(.13);
    while(B_EXT < .2){
        printf("B: %f\n", B_EXT);
        simulate(EQ_TIME, .13);
        // Measure magnetization
        results[sample_counter][0] = B_EXT;
        for(i=0; i <= NUM_L; i++)
          results[sample_counter][i+1] = calc_magnetization(i-1);

        for(cor_count = 1; cor_count < COR_TIME; cor_count++){
          simulate(1, .13);
          for(i=0; i <= NUM_L; i++)
            results[sample_counter][i+1] = calc_magnetization(i-1);
        }
        for(i=0; i <= NUM_L; i++)
          results[sample_counter][i+1] += results[sample_counter][i+1]/COR_TIME;

        sample_counter += 1;

        B_EXT += DELTA_B;
    }

    /* UNCOMMENT FOR NEGATIVE SWEEPING !!!!!!!!!!!!!!!!!!!!!!!!!!! */

    while(B_EXT > -.2){
        printf("B: %f\n", B_EXT);
        simulate(EQ_TIME, .13);
        // Measure magnetization
        results[sample_counter][0] = B_EXT;
        for(i=0; i <= NUM_L; i++)
          results[sample_counter][i+1] = calc_magnetization(i-1);

        // Take average over 1000 sweeps
        for(cor_count = 1; cor_count < COR_TIME; cor_count++){
          simulate(1, .13);
          for(i=0; i <= NUM_L; i++)
            results[sample_counter][i+1] = calc_magnetization(i-1);
        }
        for(i=0; i <= NUM_L; i++)
          results[sample_counter][i+1] += results[sample_counter][i+1]/COR_TIME;
        sample_counter += 1;
        B_EXT -= DELTA_B;
    }


    //clock_t end = clock();
    //double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("Execution time: %f \n", time_spent);

    return sample_counter;

}

double calc_TC(){
  /* Calculate the topological charge by triangulating the lattice and
     Sum over all solid angles of each three spin config and divide by 4pi.
     Solid angles computed by Berg formula.
     At site (i,j), calc solid angle for triangles [(i,j), (i-1,j), (i,j-1)],
     [(i,j), (i-1,j), (i,j+1)], [(i,j), (i,j-1), (i+1, j)], [(i,j), (i,j+1), (i+1,j)]
  */
  double solid_angle_sum = 0;
  printf("Calc tc");
  int i, j, k;
  for(i = 0; i < NUM_L; i++)
    for(j=0; j < ROWS; j++)
      for(k=0; k < COLS; k++){
        //solid_angle_sum += calc_solid_angle(lattice[i][j][k], lattice[i][(((j-1)%ROWS) + ROWS) % ROWS][k], lattice[i][j][(((k-1)%COLS) + COLS) % COLS]);
        //solid_angle_sum += calc_solid_angle(lattice[i][j][k], lattice[i][(((j-1)%ROWS) + ROWS) % ROWS][k], lattice[i][j][(k+1)%COLS]);
        //printf("Triangle (%d,%d) -> (%d,%d) -> (%d,%d): %f\n", j,k, j, (((k-1)%COLS) + COLS) % COLS, (j+1)%ROWS, k, calc_solid_angle(lattice[i][j][k], lattice[i][j][(((k-1)%COLS) + COLS) % COLS], lattice[i][(j+1)%ROWS][k]));
        solid_angle_sum += calc_solid_angle(lattice[i][j][k], lattice[i][j][(((k-1)%COLS) + COLS) % COLS], lattice[i][(j+1)%ROWS][k]);
        //solid_angle_sum += calc_solid_angle(lattice[i][j][k], lattice[i][(j+1)%ROWS][k], lattice[i][j][(k+1)%COLS]);
        //printf("Triangle (%d,%d) -> (%d,%d) -> (%d,%d): %f\n", j,k, (j)%ROWS, (k+1)%COLS, (((j-1)%ROWS) + ROWS) % ROWS, k, calc_solid_angle(lattice[i][j][k], lattice[i][(j)%ROWS][(k+1)%COLS], lattice[i][(((j-1)%ROWS) + ROWS) % ROWS][(k)%COLS]));

        solid_angle_sum += calc_solid_angle(lattice[i][j][k], lattice[i][(j)%ROWS][(k+1)%COLS], lattice[i][(((j-1)%ROWS) + ROWS) % ROWS][(k)%COLS]);

        //printf("Solid angle sum is now %f \n", solid_angle_sum);

      }

  return solid_angle_sum/(4*M_PI);

}

double calc_solid_angle(spin_t n1, spin_t n2, spin_t n3){

  gsl_complex c_temp;
  gsl_vector * n1_vec = gsl_vector_alloc(3);
  gsl_vector * n2_vec = gsl_vector_alloc(3);
  gsl_vector * n3_vec = gsl_vector_alloc(3);
  gsl_vector * n2_cross_n3 = gsl_vector_alloc(3);

  double n1_dot_n2;
  double n2_dot_n3;
  double n3_dot_n1;
  double n1_dot_n2_cross_n3;

  double rho;
  double Omega;


  gsl_vector_set(n1_vec, 0, n1.x);
  gsl_vector_set(n1_vec, 1, n1.y);
  gsl_vector_set(n1_vec, 2, n1.z);

  gsl_vector_set(n2_vec, 0, n2.x);
  gsl_vector_set(n2_vec, 1, n2.y);
  gsl_vector_set(n2_vec, 2, n2.z);

  gsl_vector_set(n3_vec, 0, n3.x);
  gsl_vector_set(n3_vec, 1, n3.y);
  gsl_vector_set(n3_vec, 2, n3.z);

  cross_product(n2_vec, n3_vec, n2_cross_n3);

  gsl_blas_ddot(n1_vec, n2_vec, &n1_dot_n2);
  gsl_blas_ddot(n2_vec, n3_vec, &n2_dot_n3);
  gsl_blas_ddot(n3_vec, n1_vec, &n3_dot_n1);
  gsl_blas_ddot(n1_vec, n2_cross_n3, &n1_dot_n2_cross_n3);

  //printf("n1_dot_n2: %f\n", n1_dot_n2);
  //printf("n2_dot_n3: %f\n", n2_dot_n3);
  //printf("n3_dot_n1: %f\n", n3_dot_n1);

  rho = sqrt(2*(1+n1_dot_n2)*(1+n2_dot_n3)*(1+n3_dot_n1));
  //printf("Rho is %f \n", rho);
  GSL_SET_COMPLEX(&c_temp, (1.0/rho)*(1 + n1_dot_n2 + n2_dot_n3 + n3_dot_n1), (1.0/rho)*n1_dot_n2_cross_n3);
  Omega = 2*GSL_IMAG(gsl_complex_log(c_temp));

  gsl_vector_free(n1_vec);
  gsl_vector_free(n2_vec);
  gsl_vector_free(n3_vec);
  gsl_vector_free(n2_cross_n3);


  return Omega;

}


void sample_mag(double* mag_vals, int sweeps){
  int n, i = 0;
  simulate(EQ_TIME, TEMP);
  for(i=0; i < sweeps; i++){
    simulate(1,TEMP);
    mag_vals[i] = calc_magnetization(-1);
    if(i%10000 == 0)
      printf("Sweep %d: %f\n", i, mag_vals[i]);
  }


    FILE *f = fopen(SIM_CONFIG, "a");
    fprintf(f, "Simulation %d (Autocorrelation)\n", SIM_NUM);
    fprintf(f, "Lattice Record %d: Size = %d\n", SIM_NUM, ROWS);
    fprintf(f,"\tD = %f\n", D);
    fprintf(f,"\tB = %f\n", B_EXT);
    fprintf(f,"\tT = %f\n", TEMP);
    fprintf(f,"\tOverrelaxation: %d\n", OVER_FLAG);
    fprintf(f,"\tSamples = %d sweeps \n", sweeps);
    fprintf(f,"\tEQ time = %d sweeps \n", EQ_TIME);

    for(n = 0; n < NUM_L; n++)
      fprintf(f, "\tJ_AF[%d] = %f\n\tK[%d] = %f\n", n, J_INTER[n], n, K[n]);
    fclose(f);

}

int mag_v_temp(double** mag_vals){
  /* Perform num_sweeps number of sweeps, measuring magnetization
  after every sweep and recording in mag_vals. */
  simulate(EQ_TIME, 5.0);
  float curr_temp = 5.0;
  int i = 0;

  mag_vals[i][0] = curr_temp;
  mag_vals[i][1] = calc_magnetization(-1);
  while(curr_temp > 0.3){
    simulate(EQ_TIME,curr_temp);
    mag_vals[i][0] = curr_temp;
    mag_vals[i++][1] = calc_magnetization(-1);
    curr_temp -= DELTA_T;

  }

  FILE *f = fopen(SIM_CONFIG, "a");
  fprintf(f, "Simulation %d: Mag v Temp\n", SIM_NUM);
  fprintf(f, "Lattice Record %d: Size = %d\n", SIM_NUM, ROWS);
  fprintf(f,"\tD = %f\n", D);
  fprintf(f,"\tB = %f\n", B_EXT);
  fprintf(f, "\tK = %f\n", K[0]);
  fprintf(f,"\tOverrelaxation: %d\n", OVER_FLAG);
  fclose(f);


  return i;

}

int suscept_v_temp(double** s_vals){
  /* Perform measuring susceptibility
  after every sweep and recording in s_vals. */
  float curr_temp = 5.0;
  int i = 0;
  int j = 0;
  double mag = 0.0;
  double mag_sum = 0.0;
  double mag_sq_sum = 0.0;


  while(curr_temp > 0.2){
    mag_sum = 0.0;
    mag_sq_sum = 0.0;
    simulate(EQ_TIME, curr_temp);
    for(j = 0; j < COR_TIME; j++){
      sweep(curr_temp);
      mag = calc_magnetization(-1);
      mag_sum += mag;
      mag_sq_sum += gsl_pow_2(mag);
    }
    s_vals[i][0] = curr_temp;
    //s_vals[i][1] = mag_sum/COR_TIME;
    //s_vals[i][2] = mag_sq_sum/COR_TIME;
    s_vals[i][1] = ((mag_sq_sum/COR_TIME) - gsl_pow_2((mag_sum/COR_TIME)))*(NUM_L*ROWS*COLS/curr_temp);
    printf("T: %f, Susceptibility: %f\n", s_vals[i][0], s_vals[i][1]);
    curr_temp -= DELTA_T;
    i++;
  }

  FILE *f = fopen(SIM_CONFIG, "a");
  fprintf(f, "Simulation %d: Suscept v Temp\n", SIM_NUM);
  fprintf(f, "Lattice Record %d: Size = %d\n", SIM_NUM, ROWS);
  fprintf(f,"\tD = %f\n", D);
  fprintf(f,"\tB = %f\n", B_EXT);
  fprintf(f, "\tK = %f\n", K[0]);

  fprintf(f,"\tOverrelaxation: %d\n", OVER_FLAG);
  fclose(f);


  return i;

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

void test_lattice_TC(){
  int i = 0;
  int j = 0;
  for(j=0; j < 5; j++){
    lattice[0][0][j].x = 0.0;
    lattice[0][0][j].y = 0.0;
    lattice[0][0][j].z = 1.0;

    lattice[0][4][j].x = 0.0;
    lattice[0][4][j].y = 0.0;
    lattice[0][4][j].z = 1.0;

  }

  for(i=1; i < 4; i++){
    lattice[0][i][0].x = 0.0;
    lattice[0][i][0].y = 0.0;
    lattice[0][i][0].z = 1.0;

    lattice[0][i][4].x = 0.0;
    lattice[0][i][4].y = 0.0;
    lattice[0][i][4].z = 1.0;

  }

  lattice[0][1][2].x = -1.0;
  lattice[0][1][2].y = 0.0;
  lattice[0][1][2].z = 0.0;

  lattice[0][3][2].x = 1.0;
  lattice[0][3][2].y = 0.0;
  lattice[0][3][2].z = 0.0;

  lattice[0][2][1].x = 0.0;
  lattice[0][2][1].y = -1.0;
  lattice[0][2][1].z = 0.0;

  lattice[0][2][3].x = 0.0;
  lattice[0][2][3].y = 1.0;
  lattice[0][2][3].z = 0.0;

  lattice[0][2][2].x = 0.0;
  lattice[0][2][2].y = 0.0;
  lattice[0][2][2].z = -1.0;

  lattice[0][1][1].x = -1.0/sqrt(2.0);
  lattice[0][1][1].y = -1.0/sqrt(2.0);
  lattice[0][1][1].z = 0.0;

  lattice[0][3][1].x = 1.0/sqrt(2.0);
  lattice[0][3][1].y = -1.0/sqrt(2.0);
  lattice[0][3][1].z = 0.0;

  lattice[0][3][3].x = 1.0/sqrt(2.0);
  lattice[0][3][3].y = 1.0/sqrt(2.0);
  lattice[0][3][3].z = 0.0;

  lattice[0][1][3].x = -1.0/sqrt(2.0);
  lattice[0][1][3].y = 1.0/sqrt(2.0);
  lattice[0][1][3].z = 0.0;

  printf("Topological charge is %f \n", calc_TC());


}
