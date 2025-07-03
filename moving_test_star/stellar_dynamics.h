#ifndef STELLAR_DYNAMICS_H
#define STELLAR_DYNAMICS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <omp.h>

// Constants
#define N 10000
#define R 10.0
#define G 1
#define SIMULATION_COUNT 10000
#define M 1.0
#define GMM 1
#define p 1.0
#define KING_W0 5.0
#define KING_RC 1.0    // Core radius
#define KING_RT 5.0    // Tidal radius (for W0=5)
#ifdef USE_KROUPA_IMF
#define M_MIN 0.08
#define M_MAX 150.0
#define ALPHA1 1.3
#define ALPHA2 2.3
#endif

enum Distribution { SPHERE = 1, PLUMMER, KING };

// Star structure
typedef struct {
    double x, y, z;
    double mass;
} Star;

// Global variables
extern Star cluster[N];
extern int option;
extern double force_array[SIMULATION_COUNT];
extern double deltaE_array[SIMULATION_COUNT];
extern double deltaLx_array[SIMULATION_COUNT];
extern double deltaLy_array[SIMULATION_COUNT];
extern double deltaLz_array[SIMULATION_COUNT];
extern double deltaL_array[SIMULATION_COUNT];
extern double F_x_array[SIMULATION_COUNT];
extern double F_y_array[SIMULATION_COUNT];
extern double F_z_array[SIMULATION_COUNT];
extern double deltaE_from_potential_array[SIMULATION_COUNT];

// Function declarations
double rand_double(gsl_rng *rng);
#ifdef USE_KROUPA_IMF
double generate_kroupa_mass(gsl_rng *rng);
#endif
void generate_test_star_position(double *x, double *y, double *z, int option, gsl_rng *rng);
void generate_star_position(Star *star, int option, gsl_rng *rng);
void generate_cluster(gsl_rng *rng);
void test_star_velocity(double x, double y, double z, double *vx, double *vy, double *vz, gsl_rng *rng);
void calculate_force(double x_test, double y_test, double z_test, double *F_x, double *F_y, double *F_z);
void calculate_angular_momentum(double x, double y, double z, double vx, double vy, double vz, 
                              double *Lx, double *Ly, double *Lz);
double calculate_energy(double x, double y, double z, double vx, double vy, double vz);
double plummer_potential(double x, double y, double z);
double sphere_potential(double x, double y, double z);
double king_potential(double x, double y, double z);
double calculate_e_wrt_potential(double x, double y, double z, double vx, double vy, double vz);
double escape_velocity(double x, double y, double z);
#endif