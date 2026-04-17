#!/bin/bash

cd flibe 

for i in $(seq 1 0.2 3.2); do # velocities, m/s

    echo "Running FLiBe Velocity Case $i"

    cd vel_${i}

    gmshToFoam tes_openfoam.msh
    checkMesh

    cd ..
    cd ..
    python3 change_patch_names.py --breeder "flibe" --velocity $i

    cd flibe/vel_${i}

    foamRun -solver incompressibleFluid

    echo "Ran FLiBe Velocity Case $i"

    cd ..

done