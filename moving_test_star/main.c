#include "stellar_dynamics.h"

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

    // Initialize OpenMPI
    int num_threads = omp_get_max_threads();
    printf("Running with %d threads\n", num_threads);
    #ifdef USE_KROUPA_IMF
    printf("Using Kroupa IMF (M_min=%.2f, M_max=%.1f)\n", M_MIN, M_MAX);
    #else
    printf("Using uniform masses (M=%.1f)\n", M);
    #endif

    gsl_rng **rngs = (gsl_rng **)malloc(num_threads * sizeof(gsl_rng *));
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        rngs[tid] = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(rngs[tid], time(NULL) + clock() + tid);
    }

    FILE *fp_all = fopen("all_data.dat", "w");
    if (fp_all == NULL) {
        fprintf(stderr, "Error opening output file\n");
        return 1;
    }
    
    double start_time = omp_get_wtime();
    
    #pragma omp parallel for
    for (int i = 0; i < SIMULATION_COUNT; i++) {
        int tid = omp_get_thread_num();
        gsl_rng *rng = rngs[tid];
        
        // Generate the cluster once per simulation
        generate_cluster(rng);
        
        double F_x = 0.0, F_y = 0.0, F_z = 0.0;
        double x_test, y_test, z_test, vx_test, vy_test, vz_test;
        double Lx_initial, Ly_initial, Lz_initial;

        generate_test_star_position(&x_test, &y_test, &z_test, option, rng);
        test_star_velocity(x_test, y_test, z_test, &vx_test, &vy_test, &vz_test, rng);
        calculate_angular_momentum(x_test, y_test, z_test, vx_test, vy_test, vz_test, &Lx_initial, &Ly_initial, &Lz_initial);
        
        calculate_force(x_test, y_test, z_test, &F_x, &F_y, &F_z);
        double E_initial = calculate_energy(x_test, y_test, z_test, vx_test, vy_test, vz_test);
        double E_from_potential_initial = calculate_e_wrt_potential(x_test, y_test, z_test, vx_test, vy_test, vz_test);
        double dt = 0.001; // dt must be small 0.01 gives unstable results.
        
        // Velocity Verlet integration
        vx_test += 0.5 * F_x / M * dt;
        vy_test += 0.5 * F_y / M * dt;
        vz_test += 0.5 * F_z / M * dt;
        
        x_test += vx_test * dt;
        y_test += vy_test * dt;
        z_test += vz_test * dt;
        
        calculate_force(x_test, y_test, z_test, &F_x, &F_y, &F_z);
        
        vx_test += 0.5 * F_x / M * dt;
        vy_test += 0.5 * F_y / M * dt;
        vz_test += 0.5 * F_z / M * dt;
        
        double Lx_final, Ly_final, Lz_final;
        calculate_angular_momentum(x_test, y_test, z_test, vx_test, vy_test, vz_test, &Lx_final, &Ly_final, &Lz_final);
        double E_final = calculate_energy(x_test, y_test, z_test, vx_test, vy_test, vz_test);
        double E_from_potential_final = calculate_e_wrt_potential(x_test, y_test, z_test, vx_test, vy_test, vz_test);
        double deltaLx = Lx_final - Lx_initial;
        double deltaLy = Ly_final - Ly_initial;
        double deltaLz = Lz_final - Lz_initial;
        double deltaL = sqrt(deltaLx*deltaLx + deltaLy*deltaLy + deltaLz*deltaLz);
        double deltaE = E_final - E_initial;
        double deltaE_from_potential = E_from_potential_final - E_from_potential_initial;

        F_x_array[i] = F_x;
        F_y_array[i] = F_y;
        F_z_array[i] = F_z;
        deltaE_array[i] = deltaE;
        deltaLx_array[i] = deltaLx;
        deltaLy_array[i] = deltaLy;
        deltaLz_array[i] = deltaLz;
        deltaL_array[i] = deltaL;
        deltaE_from_potential_array[i] = deltaE_from_potential;
    }
    
    for (int i = 0; i < SIMULATION_COUNT; i++) {
        fprintf(fp_all, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                F_x_array[i], F_y_array[i], F_z_array[i],
                deltaE_array[i], deltaL_array[i],
                deltaLx_array[i], deltaLy_array[i], deltaLz_array[i], deltaE_from_potential_array[i]);
    }
    
    fclose(fp_all);
    
    double total_time = omp_get_wtime() - start_time;
    printf("Simulation completed in %.2f seconds (%.1f sims/sec)\n", 
           total_time, SIMULATION_COUNT/total_time);
    
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        gsl_rng_free(rngs[tid]);
    }
    free(rngs);
    
    return 0;
}