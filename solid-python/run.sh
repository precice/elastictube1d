#!/bin/bash
cd ${0%/*} || exit 1    		    		# Run from this directory

# Solid participant

# Run this script in one terminal in order to start this participant.

echo "Preparing and running the Solid participant..."

python3 ./SolidSolver.py ../precice-config.xml 

echo ""