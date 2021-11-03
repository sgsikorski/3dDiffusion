# 3dDiffusion
3dDiffusion2.py

3d diffusion equation solver

Written by Scott Sikorski and Adam B. Sefkow
Last revised October 21, 2021

Still in the process of being updated, adding features, and optimizing the algorithm.

### To Execute
``` 
python3 3dDiffusion2.py
```
Process will take a while to output matplotlib plots. With current 64 cells in each dimension, 64^3 total cells, DSMC plots take about 15 minutes to compute while QDSMC plots take about 27 minutes to compute.

![Plots of 3d Diffusion Model](https://github.com/sgsikorski/3dDiffusion/endResults3dDiff.jiff)

### Features
- Creates 12 plots. Of each row, the 3 plots are slices through the middle of each axis.
- Progression of DSMC method shown from 27 particles per cell to 1000 particles per cell
- Last row is QDSMC approach with only 4 particles per cell! Super accurate for a very low memory and relatively low runtime requirements

### Why
Right now particle in cell methods are very costly in both space and runtime requirements. So we use QDSMC, an algorithm provided by Albright et al. [^1] to use less memory and improve runtime while also reducing the amount of noise that accompanies the procedure.
This method can now be used to simulate real life experiments and use the result of the simulations instead of a physical experiments. 
So combined with the runtime improvement, this method will be used for the University of Rochester Laboratory for Laser Energetics for their simulations.
And we were able to expand Albright's algorithm to 3D allowing more verstaile and realistic use.


### To Come
- Optimize the QDSMC algorithm
- User species definition (i.e, a plot of water or carbon or any molecule really)
- User defined initial conditions for density and temperature
- Euler equations used
- Velocity, temperature, pressure, and other properties computed, by use of the Euler equations

[^1] Quiet direct simulation of Eulerian fluids by Albright et al
