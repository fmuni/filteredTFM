#!/bin/sh

cd ${0%/*} || exit 1    # run from this directory

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

cleanCase

rm -r CFD/processor*
rm -r CFD/postProc*
rm -r CFD/log*
rm -r CFD/clockD*
rm -r CFD/octave/*.eps
rm -r CFD/constant/polyMesh/boundary
rm -r CFD/constant/polyMesh/faces
rm -r CFD/constant/polyMesh/neighbour
rm -r CFD/constant/polyMesh/owner
rm -r CFD/constant/polyMesh/points
