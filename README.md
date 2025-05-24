# force-distribution-on-a-star
holtsmark.c shows the force distribution on a randomly positioned test star in a cluster. Cluster under study can be run with 3 different models including homogeneous sphere, Plummer and King models. One can run the program for each of these models by changing the OPTION macro defined at the beginning of the code. Code assumes stars with equal mass and constants in the equations are normalized for simplicity.

holtsmark.c necessitates the use of GNU Scientific Library (GSL) and can be compiled with:
```
gcc holtsmark.c -lgsl -lm -o holtsmark
```

And the obtained data can be studied with plot.py

Note: plot.py ignores the outliers in the data (data points which differ significantly than other which are below 0.5% and above 99.5% in this code) as they disrupt the histogram.
