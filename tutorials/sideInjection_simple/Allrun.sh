#!/bin/sh

currDir=$PWD

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions


#Prepare mesh and BCs
cd octave
octave<generateCase.m
cd $currDir
./CreateMesh.sh


#Run the Case
runApplication decomposePar
runParallel eulerianFilteredTFM
