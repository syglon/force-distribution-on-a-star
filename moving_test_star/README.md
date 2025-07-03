# moving_test_star

This simulation models the evolution of a single test star moving through a static star cluster, using one of three potential models: **Sphere**, **Plummer**, or **King**. The goal is to study how the test star's **energy**, **angular momentum**, and **gravitational forces** vary as it interacts with the cluster.

---

## Compilation

To compile the simulation, use the Makefile provided.

- For **uniform-mass** stars:
  ```
  make
  ```
- To use the **Kroupa Initial Mass Function (IMF)** for stellar masses:
  ```
  make USE_KROUPA=1
  ```
---

## Running the Simulation

After compilation, run the simulation with one of the following cluster models:
```
  ./stellar_simulation [SPHERE | PLUMMER | KING]
```
Where:
- SPHERE: Homogeneous density sphere
- PLUMMER: Classic Plummer density profile
- KING: Tidally truncated King model

With default values for `SIMULATION_COUNT` and `N`:
- SPHERE and PLUMMER simulations take ~3 seconds
- KING model takes ~60 seconds due to more complex sampling

The output will be saved to:
  all_data.dat

---

## Analyzing Results

You can visualize the simulation output using the following Python scripts:

- To generate histograms of energy and angular momentum changes:
  ```python histograms.py```

- To analyze force distributions and compare them to theoretical models (e.g. Holtsmark or Gaussian):
  ```python forces.py```

---

In their respective folders you can see the example output histograms.
