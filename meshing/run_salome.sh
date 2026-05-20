#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

helpFunction() 
{
   echo ""
   echo "Usage: $0 -m mesh_folder -r bend_r -l length"
   echo -e "\t-m Path to mesh folder."
   echo -e "\t-r Bend radius of pipe (m)."
   echo -e "\t-l Length of pipe (m)."
   exit 1 # exit script after printing help 
}

while getopts "m:r:l:" opt
do
   case "$opt" in
      m ) mesh_folder="$OPTARG" ;;
      r ) bend_r="$OPTARG" ;;
      l ) length="$OPTARG" ;;
      ? ) helpFunction ;; # print helpFunction in case parameter is non-existent
   esac
done

# print helpFunction in case parameters are empty
if [ -z "$mesh_folder" ] || [ -z "$bend_r" ] || [ -z "$length" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# find conda base path 
CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -z "$CONDA_BASE" ]; then
    echo "ERROR: conda not found"
    exit 1
fi

source "$CONDA_BASE/etc/profile.d/conda.sh"

conda deactivate

export MESH_DIR="$(echo "$mesh_folder" | xargs)"
export BEND_R="$(echo "$bend_r" | xargs)"
export LENGTH="$(echo "$length" | xargs)"

export PYTHONWARNINGS="ignore::SyntaxWarning" # suppress salome warnings

cd /opt/salome && /opt/salome/salome -t python "$SCRIPT_DIR/one_vol_salome.py" # create one vol mesh 
cd /opt/salome && /opt/salome/salome -t python "$SCRIPT_DIR/two_vol_salome.py" # create two vol mesh 

conda activate tes-pav-env