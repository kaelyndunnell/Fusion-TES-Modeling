#!/bin/bash

helpFunction() 
{
   echo ""
   echo "Usage: $0 -p case_pathway -v v_in -b breeder -m mesh_file -r bend_radius -l length"
   echo -e "\t-p Path to pipe OpenFOAM case directory."
   echo -e "\t-b Breeder of simulation, lipb or flibe."
   echo -e "\t-v Breeder inlet velocity."
   echo -e "\t-m Pipe mesh file."
   echo -e "\t-r Bend radius."
   echo -e "\t-l Length."
   exit 1 # exit script after printing help 
}

while getopts "p:b:v:m:r:l:" opt
do
   case "$opt" in
      p ) case_pathway="$OPTARG" ;;
      b ) breeder="$OPTARG" ;;
      v ) v_in="$OPTARG" ;;
      m ) mesh_file="$OPTARG" ;;
      r ) bend_radius="$OPTARG" ;;
      l ) length="$OPTARG" ;;
      ? ) helpFunction ;; # print helpFunction in case parameter is non-existent
   esac
done

# print helpFunction in case parameters are empty
if [ -z "$case_pathway" ] || [ -z "$breeder" ] || [ -z "$v_in" ] || [ -z "$mesh_file" ] || [ -z "$bend_radius" ] || [ -z "$length" ] 
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# begin script when all parameters given 
cd $case_pathway

ideasUnvToFoam $mesh_file > log 2>&1
checkMesh > log 2>&1

cd ..
cd ..
python3 change_patch_names.py --breeder $breeder --velocity $v_in --geometry "pipe" --bend_radius $bend_radius --length $length

cd $breeder

parent_dir="$PWD"
v_in=$(echo "$v_in" | xargs)
bend_radius=$(echo "$bend_radius" | xargs)
length=$(echo "$length" | xargs)
# tank_path="$parent_dir/inlet_tank_${v_in}m_s"

cd "pipe_${v_in}m_s_r${bend_radius}_l${length}"

# mapFields $tank_path -sourceTime latestTime > log 2>&1

foamRun -solver incompressibleFluid > log 2>&1

touch tes.foam

echo "OpenFOAM pipe case completed successfully."

cd ..