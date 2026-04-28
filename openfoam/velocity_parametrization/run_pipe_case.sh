#!/bin/bash

helpFunction() 
{
   echo ""
   echo "Usage: $0 -p case_pathway -b -v v_in -t time"
   echo -e "\t-p Path to pipe OpenFOAM case directory."
   echo -e "\t-b Breeder of simulation, lipb or flibe."
   echo -e "\t-v Breeder inlet velocity."
   exit 1 # exit script after printing help 
}

while getopts "p:d:b:v:t:" opt
do
   case "$opt" in
      p ) case_pathway="$OPTARG" ;;
      b ) breeder="$OPTARG" ;;
      v ) v_in="$OPTARG" ;;
      ? ) helpFunction ;; # print helpFunction in case parameter is non-existent
   esac
done

# print helpFunction in case parameters are empty
if [ -z "$case_pathway" ] || [ -z "$breeder" ] || [ -z "$v_in" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# begin script when all parameters given 
cd $case_pathway

gmshToFoam tes_openfoam.msh > log 2>&1
checkMesh > log 2>&1

cd ..
cd ..
python3 change_patch_names.py --breeder $breeder --velocity $v_in --geometry "pipe"
echo "Wall "patch" changed to "wall""

cd $breeder

parent_dir="$PWD"
v_in=$(echo "$v_in" | xargs)
tank_path="$parent_dir/inlet_tank_${v_in}m_s"

echo $tank_path

cd "pipe_${v_in}m_s"

mapFields $tank_path -sourceTime latestTime > log 2>&1

echo done
 
foamRun -solver incompressibleFluid > log 2>&1

touch tes.foam

cd ..