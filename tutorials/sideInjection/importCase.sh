#!/bin/sh

baseCase='../sideInjection_simple'


rm -r system constant 0 octave

cp -r $baseCase/system .
cp -r $baseCase/constant .
cp -r $baseCase/0 .
cp -r $baseCase/octave .

cp blockMeshDict.orig system/blockMeshDict
cp topoSetDict.orig system/topoSetDict
cp extrudeMeshDict.orig system/extrudeMeshDict

