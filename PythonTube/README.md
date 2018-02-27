# Setup Instructions

You can use anaconda to get this example running

* install anaconda https://www.anaconda.com/download/#linux
* define ```export ANACONDA_ROOT=<path/to/anaconda/installation>```, ```export PRECICE_ROOT=<path/to/precice/folder>``` and ```export PRECICE_ROOT=<path/to/folder/containing/precice/folder>```
* run ```./precice_elastictube.sh```. This file initializes your anaconda environment ```precice_tube``` and creates the python adapter.
* activate the environment ```precice_tube``` via ```source activate precice_tube```
* go to ```cd $PRECICE_BASE/elastictube1d/PythonTube```
* run ```python FluidSolver.py precice-config.xml``` and ```python StructureSolver.py precice-config.xml``` each in one shell. Don't forget to activate the environment ```precice_tube```.
