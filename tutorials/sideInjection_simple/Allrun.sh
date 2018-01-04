#!/bin/sh

currDir=$PWD

#Prepare mesh and BCs
cd octave
octave<generateCase.m
cd $currDir
./CreateMesh.sh


#Run the Case
decomposePar
mpirun -np 4 eulerianFilteredTFM -parallel
