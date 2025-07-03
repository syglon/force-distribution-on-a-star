# force-distribution-on-a-star
holtsmark.c shows the force distribution on a randomly positioned test star in a cluster. Cluster under study can be run with 3 different models including homogeneous sphere, Plummer and King models. Code assumes stars with equal mass and constants in the equations are normalized for simplicity.

holtsmark.c necessitates the installation of GNU Scientific Library (GSL) and can be compiled with:
```
gcc holtsmark.c -lgsl -lm -o holtsmark
```
After compilation, models can be selected while executing the object file:
```
./holtsmark SPHERE # For a homogeneous sphere
./holtsmark PLUMMER # For Plummer model
```
For King model, it is highly recommended to change the SIMULATION_COUNT macro defined at the start of the code to 10000 as it may take a while to run. After the change, compile the code again and execute the object file with:
```
gcc holtsmark.c -lgsl -lm -o holtsmark
./holtsmark KING
```

And the obtained forces.dat file can be studied with plot.py
```
python plot.py
```

Notes: 
* plot.py ignores the outliers in the data (data points which differ significantly than other which are below 0.5% and above 99.5% in this code) as they disrupt the histogram.
