# Setup and Running

You can use anaconda to get this example running. Helper scripts are provided in the folder ```building```.

* install anaconda https://www.anaconda.com/download/#linux
* define ```ANACONDA_ROOT=<path/to/anaconda/installation>```, ```PRECICE_ROOT=<path/to/precice/folder>```, ```ELASTICTUBE_ROOT=<path/to/elastictube1d/folder>``` in ```precice_config.sh```.
* run ```./elastictube_install.sh```. This file initializes your anaconda environment ```precice_tube``` and creates the python adapter.
* use ```elastictube_activate.sh``` to activate the environment ```precice_tube``` via ```source elastictube_activate.sh```.
* go to ```cd $ELASTICTUBE_ROOT/PythonTube```
* run ```python FluidSolver.py precice-config.xml``` and ```python StructureSolver.py precice-config.xml``` each in one shell. Don't forget to activate the environment ```precice_tube```.
