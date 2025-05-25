#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#define SPHERE 1
#define PLUMMER 2
#define KING 3
#define OPTION PLUMMER // Available options are SPHERE, PLUMMER, KING

#define N 10000 // Star count
#define R 10.0 // Radius of the sphere
#define G 1 // Normalized gravitational constant
#define SIMULATION_COUNT 100000
#define M 1.0 // Mass of each star
#define GMM 1

double rand_double(gsl_rng *rng);
void generate_test_star_position(double *x, double *y, double *z);
void generate_star_position(double *x, double *y, double *z);

double force_array[SIMULATION_COUNT];
int count = 0;
gsl_rng *rng;

int main(){
    // GSL rng setup
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(NULL));

    FILE *fp = fopen("forces.dat", "w");
    
    for(int i =0; i < SIMULATION_COUNT; i++){
        double F_x = 0.0;
        double x_test, y_test, z_test;
        generate_test_star_position(&x_test, &y_test, &z_test);

        for (int j = 0; j < N; j++) {
            double x, y, z;
            generate_star_position(&x, &y, &z);

            double dx = (x - x_test);
            double dy = (y - y_test);
            double dz = (z - z_test);

            // Calculate the x-component of the force on a test star
            double r_squared = dx * dx + dy * dy + dz * dz;
            double r_cubed = r_squared * sqrt(r_squared);
            if (r_cubed != 0) { // Very unlikely...
                F_x += -GMM * dx / r_cubed;
            }
        }
        force_array[i] = F_x;
        count++;
        if (count % 10000 == 0) {
            printf("Processed %d simulations\n", count);
        }
    }
    for(int i = 0; i < SIMULATION_COUNT; i++) {
        fprintf(fp, "%f\n", force_array[i]);
    }
    fclose(fp);
    gsl_rng_free(rng);
    return 0;
}

double rand_double(gsl_rng *rng){
    return gsl_rng_uniform(rng);
}

// Generate test star position according to OPTION
void generate_test_star_position(double *x, double *y, double *z) {
    #if OPTION == SPHERE
        double x_test, y_test, z_test;
        do {
            x_test = 2*rand_double(rng) - 1;
            y_test = 2*rand_double(rng) - 1;
            z_test = 2*rand_double(rng) - 1;
        } while (x_test*x_test + y_test*y_test + z_test*z_test > 1.0);
    #elif OPTION == PLUMMER
        double a_test = 1.0;
        double X = rand_double(rng);
        double r_test = a_test / sqrt(pow(X, -2.0/3.0) - 1.0);
        double theta_test = acos(1 - 2 * rand_double(rng));
        double phi_test = 2 * M_PI * rand_double(rng);
        *x = r_test * sin(theta_test) * cos(phi_test);
        *y = r_test * sin(theta_test) * sin(phi_test);
        *z = r_test * cos(theta_test);
    #elif OPTION == KING
        double W0 = 5.0; // King model parameter
        double r_c = 1.0; // Core radius
        double r_t = R;   // Tidal radius
        
        double sigma = sqrt(G * M * N / r_t); // Velocity dispersion
        double rho0 = 1.0;
        
        // Rejection sampling for test star position
        double r_test, psi, rho_ratio, y_dir;
        while (1) {
            // Trial radius
            r_test = r_t * pow(rand_double(rng), 1.0/3.0);
            psi = W0 * sigma*sigma * (1.0/sqrt(1.0 + (r_test/r_c)*(r_test/r_c)) - 1.0/sqrt(1.0 + (r_t/r_c)*(r_t/r_c)));
            
            // Density
            rho_ratio = pow(1.0 + (r_test/r_c)*(r_test/r_c), -1.5) * exp(psi/(sigma*sigma));
            
            // Accept/reject
            y_dir = rand_double(rng);
            if (y_dir < rho_ratio) break;
        }
        // Random direction for test star
        double theta_test = acos(1 - 2 * rand_double(rng));
        double phi_test = 2 * M_PI * rand_double(rng);
        
        *x = r_test * sin(theta_test) * cos(phi_test);
        *y = r_test * sin(theta_test) * sin(phi_test);
        *z = r_test * cos(theta_test);
    #endif    
 
}

// Generate cluster star position according to OPTION
void generate_star_position(double *x, double *y, double *z) {
    #if OPTION == SPHERE
        double r = R * pow(rand_double(rng), 1.0/3.0);
    #elif OPTION == PLUMMER
        double a = 1.0;
        double X = rand_double(rng);
        double r = a / sqrt(pow(X, -2.0/3.0) - 1.0);
    #elif OPTION == KING
        // King model parameters
        double W0 = 5.0;
        double r_c = 1.0;    // Core radius
        double r_t = R;      // Tidal radius
        
        double sigma = sqrt(G * M * N / r_t); // Velocity dispersion
        double rho0 = 1.0;
        
        // Rejection sampling
        double r, psi, rho_ratio, y_dir;
        while (1) {
            // Trial radius
            r = r_t * pow(rand_double(rng), 1.0/3.0);
            psi = W0 * sigma*sigma * (1.0/sqrt(1.0 + (r/r_c)*(r/r_c)) - 1.0/sqrt(1.0 + (r_t/r_c)*(r_t/r_c)));
            
            // Density
            rho_ratio = pow(1.0 + (r/r_c)*(r/r_c), -1.5) * exp(psi/(sigma*sigma));
            
            // Accept/reject
            y_dir = rand_double(rng);
            if (y_dir < rho_ratio) break;
        }
    #endif

    double theta = acos(1 - 2 * rand_double(rng));
    double phi = 2 * M_PI * rand_double(rng);

    *x = r * sin(theta) * cos(phi);
    *y = r * sin(theta) * sin(phi);
    *z = r * cos(theta);
}
