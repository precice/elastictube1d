#!/bin/bash
cd ${0%/*} || exit 1    		    		# Run from this directory

# Fluid participant

# Run this script in one terminal in order to start this participant.

echo "Preparing and running the Fluid participant..."

python3 ./FluidSolver.py ../precice-config.xml 

echo ""