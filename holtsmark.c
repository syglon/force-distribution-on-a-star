#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>

#define N 10000           // Star count
#define R 10.0            // Radius of the sphere
#define G 1               // Normalized gravitational constant
#define SIMULATION_COUNT 100000 // Number of simulations
#define M 1.0             // Mass of each star
#define GMM 1

enum Distribution { SPHERE = 1, PLUMMER, KING };
int option = SPHERE;

double rand_double(gsl_rng *rng);
void generate_test_star_position(double *x, double *y, double *z, int option);
void generate_star_position(double *x, double *y, double *z, int option);

double force_array[SIMULATION_COUNT];
int count = 0;
gsl_rng *rng;

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s [SPHERE|PLUMMER|KING]\n", argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "SPHERE") == 0) option = SPHERE;
    else if (strcmp(argv[1], "PLUMMER") == 0) option = PLUMMER;
    else if (strcmp(argv[1], "KING") == 0) option = KING;
    else {
        fprintf(stderr, "Invalid option: %s\nValid options are: SPHERE, PLUMMER, KING\n", argv[1]);
        return 1;
    }

    // GSL rng setup
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(NULL));

    FILE *fp = fopen("forces.dat", "w");
    
    for (int i = 0; i < SIMULATION_COUNT; i++) {
        double F_x = 0.0;
        double x_test, y_test, z_test;
        generate_test_star_position(&x_test, &y_test, &z_test, option);

        for (int j = 0; j < N; j++) {
            double x, y, z;
            generate_star_position(&x, &y, &z, option);

            double dx = x - x_test;
            double dy = y - y_test;
            double dz = z - z_test;

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

    for (int i = 0; i < SIMULATION_COUNT; i++) {
        fprintf(fp, "%f\n", force_array[i]);
    }
    fclose(fp);
    gsl_rng_free(rng);
    return 0;
}

double rand_double(gsl_rng *rng) {
    return gsl_rng_uniform(rng);
}

// Generate test star position according to OPTION
void generate_test_star_position(double *x, double *y, double *z, int option) {
    if (option == SPHERE) {
        do {
            *x = 2 * rand_double(rng) - 1;
            *y = 2 * rand_double(rng) - 1;
            *z = 2 * rand_double(rng) - 1;
        } while ((*x) * (*x) + (*y) * (*y) + (*z) * (*z) > 1.0);
    } else if (option == PLUMMER) {
        double a_test = 1.0;
        double X = rand_double(rng);
        double r_test = a_test / sqrt(pow(X, -2.0 / 3.0) - 1.0);
        double theta_test = acos(1 - 2 * rand_double(rng));
        double phi_test = 2 * M_PI * rand_double(rng);
        *x = r_test * sin(theta_test) * cos(phi_test);
        *y = r_test * sin(theta_test) * sin(phi_test);
        *z = r_test * cos(theta_test);
    } else if (option == KING) {
        double W0 = 5.0;
        double r_c = 1.0;
        double r_t = R;
        double sigma = sqrt(G * M * N / r_t);

        double r_test, psi, rho_ratio, y_dir;
        while (1) {
            r_test = r_t * pow(rand_double(rng), 1.0 / 3.0);
            psi = W0 * sigma * sigma * (1.0 / sqrt(1.0 + (r_test / r_c) * (r_test / r_c)) - 1.0 / sqrt(1.0 + (r_t / r_c) * (r_t / r_c)));
            rho_ratio = pow(1.0 + (r_test / r_c) * (r_test / r_c), -1.5) * exp(psi / (sigma * sigma));
            y_dir = rand_double(rng);
            if (y_dir < rho_ratio) break;
        }

        double theta_test = acos(1 - 2 * rand_double(rng));
        double phi_test = 2 * M_PI * rand_double(rng);

        *x = r_test * sin(theta_test) * cos(phi_test);
        *y = r_test * sin(theta_test) * sin(phi_test);
        *z = r_test * cos(theta_test);
    }
}

// This function differ from the previous one only for the SPHERE option. Stars are uniformly distributed in a sphere.
void generate_star_position(double *x, double *y, double *z, int option) {
    double r = 0.0;
    if (option == SPHERE) {
        r = R * pow(rand_double(rng), 1.0 / 3.0);
    } else if (option == PLUMMER) {
        double a = 1.0;
        double X = rand_double(rng);
        r = a / sqrt(pow(X, -2.0 / 3.0) - 1.0);
    } else if (option == KING) {
        double W0 = 5.0;
        double r_c = 1.0;
        double r_t = R;
        double sigma = sqrt(G * M * N / r_t);

        double psi, rho_ratio, y_dir;
        while (1) {
            r = r_t * pow(rand_double(rng), 1.0 / 3.0);
            psi = W0 * sigma * sigma * (1.0 / sqrt(1.0 + (r / r_c) * (r / r_c)) - 1.0 / sqrt(1.0 + (r_t / r_c) * (r_t / r_c)));
            rho_ratio = pow(1.0 + (r / r_c) * (r / r_c), -1.5) * exp(psi / (sigma * sigma));
            y_dir = rand_double(rng);
            if (y_dir < rho_ratio) break;
        }
    }

    double theta = acos(1 - 2 * rand_double(rng));
    double phi = 2 * M_PI * rand_double(rng);
    *x = r * sin(theta) * cos(phi);
    *y = r * sin(theta) * sin(phi);
    *z = r * cos(theta);
}
