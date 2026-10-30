Metropolis MC of a surface step

Monte Carlo simulation of the growth of a step on a Kossel (100) surface (SOS model), with Metropolis, n-fold way and particle-number-changing dynamics.

Files

- simulation_step.f90: Fortran module with the Monte Carlo engine (moves, Metropolis variants, n-fold way, autocorrelation, block averaging) compiled on the fly via f2py_jit.
- mc_edge.ipynb: runs the simulation for a chosen (L, T, N, move), plots the raw time evolution and saves the results to a .npz file.
- analisi_mc.ipynb: loads one or more saved runs and produces the analysis plots.

REQUIREMENTS:

Python, a Fortran compiler (gfortran) and NumPy/Matplotlib (in requirements).

sudo apt install gfortran                                      -> installation of gfortran

f2py_jit to compile the Fortran module from the notebook.

Create a virtual environment:

python -m venv env (env is the name of the virtual enviroment) -> creation of the enviroment
. env/bin/activate                                             -> activation of the enviroment
pip install -r requirements.txt                                -> installation of needed packages to run the simulations


Run mc_edge.ipynb first to generate data, then analisi_mc.ipynb to analyze it.
