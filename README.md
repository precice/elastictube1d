---
title: 1D Elastic Tube
permalink: unknown
keywords: OpenFOAM, python
summary: The 1D Elastic Tube is a FSI case, that consists of an internal flow in a flexible tube. The flow is unsteady and incompressible. This tutorial can be run using C++ or python.  Running the simulation takes a few minutes. 
---


## Setup

An accurate descrption of the setup is described [here] (https://github.com/precice/precice/wiki/1D-elastic-tube:-Case-Description). 

## Available solvers

Both fluid and solid participant are supported in:

* *cxx*
* *python*: This version is realized using the Python API for preCICE. Check [this entry in the preCICE wiki](https://github.com/precice/precice/wiki/1D-elastic-tube-using-the-Python-API) for more information on this example.


## Running the Simulation 

Before you start: A working known combination for the input parameters is N=100, tau = 0.01, kappa = 100. Other parameters like the simulation's end time, you can modify them in the precice-config.xml.

### C++

Open two separate terminals and start each participant by calling the respective run script. 

Serial run:

```
cd fluid-cpp
./FluidSolver precice-config.xml N tau kappa
```
and
```
cd solid-cpp
./SolidSolver precice-config.xml N
```
 
Parallel run:

```
cd fluid-cpp
mpiexec -np <nproc> ./FluidSolver precice-config.xml N tau kappa -parallel
```
and
```
cd solid-cpp
mpiexec -np <nproc> ./SolidSolver precice-config.xml N -parallel
```

Note that you first need to build the scripts `FluidSolver` and `SolidSolver`. To do so, run this in the root folder:

```
cmake .
make all
```

### python

Open two separate terminals and start each participant by calling the respective run script. Only serial run:

```
cd fluid-python
./FluidSolver precice-config.xml N tau kappa
```
and
```
cd solid-python
./SolidSolver precice-config.xml N
```

**Optional:** Visualization and video output of the fluid participant can be triggered via the options `--enable-plot` and `--write-video` of `FluidSolver.py`. 

## Post-processing

The visualization of results using python can be selected as an option during execution. 

In case of running the C++ version, you can visualize the results with the `Postproc/fluid.py` script:

```bash
$ python3 Postproc/fluid.py <quantity> Postproc/<prefix>
```
Note the required arguments specifying which quantity to plot (`pressure`, `velocity` or `diameter`) and a name prefix for the target vtk files.
For example, to plot the diameter using the default prefix for vtk files, we execute:
```bash
$ python3 Postproc/fluid.py diameter Postproc/out_fluid_
```
![FSI3 setup](images/diameter.png)

If you are run the case in parallel, you can visualize the results caclulated by one rank (eg. rank 0) as follows:

```bash
$ python3 Postproc/fluid.py diameter Postproc/out_fluid0_
```

An image of this diameter plot can be found in the `/images` folder.


## References

[1] M. Mehl, B. Uekermann, H. Bijl, D. Blom, B. Gatzhammer, and A. van Zuijlen.
Parallel coupling numerics for partitioned fluid-structure interaction simulations. CAMWA, 2016.  
[2] J. Degroote, P. Bruggeman, R. Haelterman, and J. Vierendeels. Stability of a coupling technique
for partitioned solvers in FSI applications. Computers & Structures, 2008.




