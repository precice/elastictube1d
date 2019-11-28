# preCICE example for FSI: 1D elastic tube

## First Step: Choose a version

This tutorial comes in two distinct versions, one written in C++ and one in Python. The code for each are located in the folders `cxx` and `python` respectively.

Choose a version, then navigate down to the corresponding sections for [Python](#python-version) and [C++](#c-version) instructions.

---
## Python version

This version is realized using the Python API for preCICE. Check [this entry in the preCICE wiki](https://github.com/precice/precice/wiki/1D-elastic-tube-using-the-Python-API) for more information on this example.

### Requirements

- [preCICE](https://github.com/precice/precice/wiki/Get-preCICE)

- [preCICE Python Bindings](https://github.com/precice/precice/blob/develop/src/precice/bindings/python/README.md).

**Note:** these requirements are specific to the Python variant of the code. If you wish to run the C++ version, you require different packages (see [this section](#c-version)).

### How to run

0. Clone this repository:
    ```bash
    $ git clone https://github.com/precice/elastictube1d.git
    ```

1. Navigate into the `python` directory and execute the `Allrun` script to run both solver participants directly:
    ```bash
    $ cd python/ && ./Allrun
    ```

2. After the script exits, you can view the output `.vtk` files in the `VTK` folder located in the root directory.

3. To clean up the log files and vtk files created during a run, execute the `Allclean` script.
    ```bash
    $ ./Allclean
    ```

---
## C++ version

Check [this preCICE wiki page](https://github.com/precice/precice/wiki/Example-for-FSI:-1D-elastic-tube) for a detailed description of this tutorial. For more information see [1]. Elastictube scenario taken from [2].

### Requirements

- [preCICE](https://github.com/precice/precice/wiki/Get-preCICE)

- [LAPACK](http://performance.netlib.org/lapack/#_lapack_version_3_8_0_2). On Ubuntu-like Linux distributions, you can also install via executing:
  ```bash
  $ sudo apt-get install liblapack-dev
  ```

### How to run

0. Clone this repository:
    ```bash
    $ git clone https://github.com/precice/elastictube1d.git
    ```

1. Navigate into the `cxx` directory and build the Makefiles:
  ```bash
  $ cd cxx/ && cmake .
  ```

2. Make the tutorial:
  ```bash
  $ make all
  ```

3. After successful compilation, you can now launch preset configuration by calling the `Allrun` script located in the current folder:
  ```bash
  $ ./Allrun
  ```
  Results will be stored as .vtk files in the `cxx/Postproc` folder.

  **Optional:** You can visualize the results with the `Postproc/fluid.py` script:
  ```bash
  $ python Postproc/fluid.py <quantity> Postproc/<prefix>
  ```
  Note the required arguments specifying which quantity to plot (`pressure`, `velocity` or `diameter`) and a name prefix for the target vtk files.
  For example, to plot the diameter using the default prefix for vtk files, we execute:
  ```bash
  $ python Postproc/fluid.py diameter Postproc/out_fluid_
  ```
  An image of this diameter plot can be found in the `cxx/example` folder.

  **Alternative:**: If you wish to run the parallel versions of each solver, run the `Allrun_parallel` script instead. Note that no vtk output is generated for this solver configuration!

5. To quickly clean the folder of log files and results from previous runs, execute `Allclean`:
  ```bash
  $ ./Allclean
  ```


**Note:** The tutorial can also be run manually by launching both participants by hand. See [this preCICE wiki page](https://github.com/precice/precice/wiki/Running-the-1D-elastic-tube-example) for instructions.

---
## References

[1] M. Mehl, B. Uekermann, H. Bijl, D. Blom, B. Gatzhammer, and A. van Zuijlen.
Parallel coupling numerics for partitioned fluid-structure interaction simulations. CAMWA, 2016.  
[2] J. Degroote, P. Bruggeman, R. Haelterman, and J. Vierendeels. Stability of a coupling technique
for partitioned solvers in FSI applications. Computers & Structures, 2008.
