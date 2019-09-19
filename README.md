# preCICE example for FSI: 1D elastic tube

## Requirements

- [preCICE](https://github.com/precice/precice/wiki/Get-preCICE)

- [preCICE Python Bindings](https://github.com/precice/precice/blob/develop/src/precice/bindings/python/README.md).

**Note:** these requirements are specific to the **Python variant** of the code. If you wish to run the C++ version, you require different packages (see [this section](#c++-variant)).

## Quick Start

1. Clone this repository:
    ```bash
    $ git clone https://github.com/Eder-K/elastictube1d.git
    ```
- Execute the `Allrun` script to directly run both solver participants directly:
    ```bash
    $ ./Allrun
    ```
    The solvers will then run a preset configuration without any further input.

- After the script exits, you can view the output `.vtk` files in the `VTK` folder located in the root directory.

- To clean up the log files and vtk files created during a run, execute the `Allclean` script.
    ```bash
    $ ./Allclean
    ```


## Contents

Check [this preCICE wiki page](https://github.com/precice/precice/wiki/Example-for-FSI:-1D-elastic-tube) for a detailed description of this tutorial. For more information see [1]. Elastictube scenario taken from [2].

This repository contains two variants of the elastictube1d tutorial, written in C++ and Python. These are located in the subfolders `elastictube1d-cxx` and `elastictube1d-python` respectively.

The script files `Allclean` and `Allrun` in the root directory provide a shortcut to run the Python version.
Similar scripts `Allclean-cxx`, `Allrun-cxx`, `Allrun_parallel-cxx` can be found in the `elastictube1d-cxx` folder and serve as shortcuts to run the C++ tutorial, both in serial and parallel modes.

### C++ variant

#### Requirements (C++)

- [preCICE](https://github.com/precice/precice/wiki/Get-preCICE)

- [LAPACK](http://performance.netlib.org/lapack/#_lapack_version_3_8_0_2). On Ubuntu-like Linux distributions, you can also install via executing:
  ```bash
  $ sudo apt-get install liblapack-dev
  ```

#### How to execute (C++)

1. Navigate into the `elastictube1d-cxx` folder and build the Makefiles:
```bash
$ cmake .
```
2. Make the tutorial:
```bash
$ make all
```
3. After successful compilation, you can now launch preset configuration by calling the `Allrun-cxx` script located in the current folder:
```bash
$ ./Allrun-cxx
```
Results will be stored in the `elastictube1d-cxx/Postproc` folder.

**Alternative**: If you wish to run the parallel versions of each solver, run the `Allrun_parallel-cxx` script instead.


4. To quickly clean the folder of log files and results from previous runs, execute `Allclean-cxx`:
```bash
$ ./Allclean
```

The tutorial can also be run manually by launching both participants by hand. See [this preCICE wiki page](https://github.com/precice/precice/wiki/Running-the-1D-elastic-tube-example) for instructions.


### Python variant

To run, simply follow the [Requirements](#requirements) and [Quick Start](#quick-start) from the top.

This version is realized using the Python API for preCICE. Check [this entry in the preCICE wiki](https://github.com/precice/precice/wiki/1D-elastic-tube-using-the-Python-API) for more information on this example.


## References

[1] M. Mehl, B. Uekermann, H. Bijl, D. Blom, B. Gatzhammer, and A. van Zuijlen.
Parallel coupling numerics for partitioned fluid-structure interaction simulations. CAMWA, 2016.  
[2] J. Degroote, P. Bruggeman, R. Haelterman, and J. Vierendeels. Stability of a coupling technique
for partitioned solvers in FSI applications. Computers & Structures, 2008.
