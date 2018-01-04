#!/bin/sh
./Allclean

blockMesh
topoSet
createPatch -overwrite
extrudeMesh
