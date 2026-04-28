#!/bin/bash

helpFunction() 
{
   echo ""
   echo "Usage: $0 -p case_pathway -b -v v_in"
   echo -e "\t-p Path to OpenFOAM case directory."
   echo -e "\t-b Breeder of simulation, lipb or flibe."
   echo -e "\t-v Breeder inlet velocity."
   exit 1 # exit script after printing help 
}

while getopts "p:b:v:" opt
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

gmshToFoam inlet_tank.msh > log 2>&1
checkMesh > log 2>&1

cd ..
cd ..
python3 change_patch_names.py --breeder $breeder --velocity $v_in --geometry "inlet_tank"

cd $breeder

v_in=$(echo "$v_in" | xargs)
cd "inlet_tank_${v_in}m_s"

foamRun -solver incompressibleFluid > log 2>&1

touch tes.foam

echo "OpenFOAM tank case completed successfully."

cd ..