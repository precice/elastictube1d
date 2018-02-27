# 1d python tube
cd $PRECICE_BASE
git clone https://github.com/precice/elastictube1d/PythonTube elastictube1d

# extend precice environment with necessary python modules
conda env create --force -f elastictube1d/PythonTube/precice_tube.yml

source $ANACONDA_ROOT/bin/activate precice_tube

# build python interface
cd $PRECICE_ROOT/src/precice/adapters/python
python setup.py build_ext --inplace || exit
