#!/bin/sh

currDir=$PWD

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

#Run the Case
runApplication blockMesh
runApplication setFields
runApplication decomposePar
runParallel eulerianFilteredTFM
