#include "stellar_dynamics.h"

// Global variables definitions
Star cluster[N];
int option = SPHERE;
double force_array[SIMULATION_COUNT];
double deltaE_array[SIMULATION_COUNT];
double deltaLx_array[SIMULATION_COUNT];
double deltaLy_array[SIMULATION_COUNT];
double deltaLz_array[SIMULATION_COUNT];
double deltaL_array[SIMULATION_COUNT];
double F_x_array[SIMULATION_COUNT];
double F_y_array[SIMULATION_COUNT];
double F_z_array[SIMULATION_COUNT];
double deltaE_from_potential_array[SIMULATION_COUNT];
double M_cluster = 0.0;

#ifdef USE_KROUPA_IMF
// Kroupa IMF parameters
#define M_MIN 0.08    // Minimum stellar mass (Solar masses)
#define M_MAX 150.0   // Maximum stellar mass (Solar masses)
#define ALPHA1 1.3    // Exponent for m < 0.5 Msun
#define ALPHA2 2.3    // Exponent for m >= 0.5 Msun

double generate_kroupa_mass(gsl_rng *rng) {
    double m;
    double xi1 = pow(0.5, -ALPHA1 + 1) - pow(M_MIN, -ALPHA1 + 1);
    double xi2 = pow(M_MAX, -ALPHA2 + 1) - pow(0.5, -ALPHA2 + 1);
    double xi_total = xi1 + xi2;
    
    double u = rand_double(rng) * xi_total;
    
    if (u <= xi1) {
        m = pow(u * (-ALPHA1 + 1) + pow(M_MIN, -ALPHA1 + 1), 1.0/(-ALPHA1 + 1));
    } else {
        m = pow((u - xi1) * (-ALPHA2 + 1) + pow(0.5, -ALPHA2 + 1), 1.0/(-ALPHA2 + 1));
    }
    
    return m;
}
#endif

// Calculate Plummer potential at position (x,y,z)
double plummer_potential(double x, double y, double z) {
    double r = sqrt(x*x + y*y + z*z);
    return -G * M_cluster / sqrt(r*r + p*p);
}

// Calculate Sphere potential at position (x,y,z)
double sphere_potential(double x, double y, double z) {
    double r = sqrt(x*x + y*y + z*z);
    if (r <= R) {
        return -0.5 * G * M_cluster * (3.0 - (r*r)/(R*R)) / R;
    } else {
        return -G * M_cluster / r;
    }
}

double king_potential(double x, double y, double z) {
    double sigma = sqrt(G * M_cluster / KING_RT);
    double r = sqrt(x*x + y*y + z*z);
    double b = r / KING_RC;

    double phi0 = G * M_cluster / sqrt(KING_RT*KING_RT + KING_RC*KING_RC);
    double phi_r = phi0 - G * M_cluster / sqrt(r*r + KING_RC*KING_RC);
    
    return (r <= KING_RT) ? phi_r : 0.0;
}

// Calculate energy using mean-field potential
double calculate_e_wrt_potential(double x, double y, double z, double vx, double vy, double vz) {
    double kinetic = 0.5 * M * (vx*vx + vy*vy + vz*vz);
    double potential;

    switch (option) {
        case SPHERE:
            potential = M * sphere_potential(x, y, z);
            break;
        case PLUMMER:
            potential = M * plummer_potential(x, y, z);
            break;
        case KING:
            potential = M * king_potential(x, y, z);
            break;
        default:
            potential = 0.0;
    }
    return kinetic + potential;
}

double calculate_energy(double x, double y, double z, double vx, double vy, double vz) {
    double kinetic = 0.5 * M * (vx*vx + vy*vy + vz*vz);
    double potential = 0.0;

    for (int i = 0; i < N; i++) {
        double dx = cluster[i].x - x;
        double dy = cluster[i].y - y;
        double dz = cluster[i].z - z;
        double r = sqrt(dx*dx + dy*dy + dz*dz);

        if (r > 0) {
            potential += -G * cluster[i].mass * M / r;
        }
    }
    return kinetic + potential;
}

void generate_cluster(gsl_rng *rng) {
    double total_mass = 0.0;

    #pragma omp parallel for reduction(+:total_mass)
    for (int i = 0; i < N; i++) {
        generate_star_position(&cluster[i], option, rng);
        total_mass += cluster[i].mass;
    }

    M_cluster = total_mass;
}

double rand_double(gsl_rng *rng) {
    return gsl_rng_uniform(rng);
}

void generate_test_star_position(double *x, double *y, double *z, int option, gsl_rng *rng) {
    if (option == SPHERE) {
        do {
            *x = R * (2 * rand_double(rng) - 1);
            *y = R * (2 * rand_double(rng) - 1);
            *z = R * (2 * rand_double(rng) - 1);
        } while ((*x)*(*x) + (*y)*(*y) + (*z)*(*z) > R*R);
    } 
    else if (option == PLUMMER) {
        double a = 1.0;
        double X = rand_double(rng);
        double r = a / sqrt(pow(X, -2.0/3.0) - 1.0);
        double theta = acos(1 - 2 * rand_double(rng));
        double phi = 2 * M_PI * rand_double(rng);
        *x = r * sin(theta) * cos(phi);
        *y = r * sin(theta) * sin(phi);
        *z = r * cos(theta);
    }
    else if (option == KING) {
        double u;
        do {
            u = rand_double(rng);
        } while (u >= 1.0); // Avoid division by 0

        double r = KING_RC * sqrt(pow(1.0 - u, -2.0/3.0) - 1.0);
        if (r > KING_RT) r = KING_RT;

        double theta = acos(1 - 2 * rand_double(rng));
        double phi = 2 * M_PI * rand_double(rng);

        *x = r * sin(theta) * cos(phi);
        *y = r * sin(theta) * sin(phi);
        *z = r * cos(theta);
    }
}

void generate_star_position(Star *star, int option, gsl_rng *rng) {
    double r = 0.0;
    
    if (option == SPHERE) {
        r = R * pow(rand_double(rng), 1.0/3.0);
    } 
    else if (option == PLUMMER) {
        double a = 1.0;
        double X = rand_double(rng);
        r = a / sqrt(pow(X, -2.0/3.0) - 1.0);
    }
    else if (option == KING) {
        double W0 = 5.0;
        double r_c = 1.0;
        double r_t = R;
        double sigma = sqrt(G * M * N / r_t);
        
        int max_iter = 100000;
        int iter = 0;
        
        while (iter++ < max_iter) {
            r = r_t * pow(rand_double(rng), 1.0/3.0);
            double rc_term = sqrt(1.0 + (r/r_c)*(r/r_c));
            double rt_term = sqrt(1.0 + (r_t/r_c)*(r_t/r_c));
            
            double psi = W0 * sigma * sigma * (1.0/rc_term - 1.0/rt_term);
            double rho_ratio = pow(1.0 + (r/r_c)*(r/r_c), -1.5) * exp(psi/(sigma*sigma));
            
            if (rand_double(rng) < rho_ratio) break;
            
            if (iter == max_iter) {
                r = r_t * pow(rand_double(rng), 1.0/3.0);
                break;
            }
        }
    }

    double theta = acos(1 - 2 * rand_double(rng));
    double phi = 2 * M_PI * rand_double(rng);
    
    star->x = r * sin(theta) * cos(phi);
    star->y = r * sin(theta) * sin(phi);
    star->z = r * cos(theta);
    
    #ifdef USE_KROUPA_IMF
        star->mass = generate_kroupa_mass(rng);
    #else
        star->mass = M;
    #endif
}

void test_star_velocity(double x, double y, double z, double *vx, double *vy, double *vz, gsl_rng *rng) {
    if (option == SPHERE) {
        double escape_vel = escape_velocity(x, y, z);
        double v = escape_vel * pow(rand_double(rng), 1.0/3.0);
        double theta = acos(1 - 2 * rand_double(rng));
        double phi = 2 * M_PI * rand_double(rng);
        *vx = v * sin(theta) * cos(phi);
        *vy = v * sin(theta) * sin(phi);
        *vz = v * cos(theta);
    }
    else if (option == PLUMMER) {
        double v_esc = escape_velocity(x, y, z);
        
        double v, f;
        do {
            v = v_esc * rand_double(rng);
            f = v*v * pow(1.0 - (v*v)/(v_esc*v_esc), 3.5);
        } while (rand_double(rng) > f);
        
        double theta = acos(1 - 2 * rand_double(rng));
        double phi = 2 * M_PI * rand_double(rng);
        *vx = v * sin(theta) * cos(phi);
        *vy = v * sin(theta) * sin(phi);
        *vz = v * cos(theta);
    }
    else if (option == KING) {
        double sigma = sqrt(G * M_cluster / KING_RT);
        double r = sqrt(x*x + y*y + z*z);
        double psi = -king_potential(x, y, z);
        double v_esc = sqrt(2.0 * psi);

        double v, f, f_max, E_over_sigma2;
        int max_iter = 100000;
        int iter = 0;

        f_max = exp(psi / (sigma * sigma)) - 1.0;

        do {
            iter++;
            if (iter > max_iter) {
                fprintf(stderr, "Warning: King velocity sampling failed at r=%.3f after %d tries\n", r, max_iter);
                v = sigma * sqrt(-2.0 * log(rand_double(rng)));
                break;
            }

            v = v_esc * rand_double(rng);
            double E = psi - 0.5 * v * v;
            E_over_sigma2 = E / (sigma * sigma);

            if (E_over_sigma2 < 0.0) continue;

            f = exp(E_over_sigma2) - 1.0;
        } while (rand_double(rng) > f / f_max);

        // Random isotropic direction
        double theta_v = acos(1 - 2 * rand_double(rng));
        double phi_v = 2 * M_PI * rand_double(rng);

        *vx = v * sin(theta_v) * cos(phi_v);
        *vy = v * sin(theta_v) * sin(phi_v);
        *vz = v * cos(theta_v);
    }

}

void calculate_force(double x_test, double y_test, double z_test, double *F_x, double *F_y, double *F_z) {
    *F_x = 0.0;
    *F_y = 0.0;
    *F_z = 0.0;
    
    for (int i = 0; i < N; i++) {
        double dx = cluster[i].x - x_test;
        double dy = cluster[i].y - y_test;
        double dz = cluster[i].z - z_test;
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        
        if (r > 0) {
            double F = -G * cluster[i].mass * M / (r*r);
            *F_x += F * dx/r;
            *F_y += F * dy/r;
            *F_z += F * dz/r;
        }
    }
}

void calculate_angular_momentum(double x, double y, double z, double vx, double vy, double vz,
                              double *Lx, double *Ly, double *Lz) {
    *Lx = M * (y * vz - z * vy);
    *Ly = M * (z * vx - x * vz);
    *Lz = M * (x * vy - y * vx);
}

double escape_velocity(double x, double y, double z) {
    double phi = 0.0;
    switch (option) {
        case SPHERE:
            phi = sphere_potential(x, y, z);
            break;
        case PLUMMER:
            phi = plummer_potential(x, y, z);
            break;
        case KING:
            phi = king_potential(x, y, z);
            break;
        default:
            phi = 0.0;
    }
    return sqrt(-2.0 * phi);
}