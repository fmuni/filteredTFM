#!/bin/sh

currDir=$PWD

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

./Allclean.sh

#Run the Case
runApplication blockMesh
#runApplication decomposePar
#runParallel eulerianFilteredTFM
runApplication eulerianFilteredTFM
