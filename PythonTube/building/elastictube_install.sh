#!/bin/bash
source ./precice_config.sh
export ELASTICTUBE_ROOT

# create a conda environment with necessary modules
conda env create --force -f $ELASTICTUBE_ROOT/PythonTube/building/precice_tube.yml

# activate conda environment and set environmental variables
source elastictube_activate.sh

# build python interface
cd $PRECICE_ROOT/src/precice/adapters/python
python setup.py build_ext --inplace || exit

cd $ELASTICTUBE_ROOT/PythonTube/building
source deactivate
