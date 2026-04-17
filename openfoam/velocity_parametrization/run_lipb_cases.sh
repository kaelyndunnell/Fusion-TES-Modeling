#!/bin/bash

cd lipb 

for i in $(seq 0.1 0.1 1); do # velocities, m/s

    echo "Running LiPb Velocity Case $i"

    cd vel_${i}

    gmshToFoam tes_openfoam.msh
    checkMesh

    cd ..
    cd ..
    python3 change_patch_names.py --breeder "lipb" --velocity $i
    echo "Wall "patch" changed to "wall""

    cd lipb/vel_${i}

    foamRun -solver incompressibleFluid

    echo "Ran LiPb Velocity Case $i"

    cd ..

done