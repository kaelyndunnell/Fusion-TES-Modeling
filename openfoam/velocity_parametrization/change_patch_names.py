import os
import shutil
import argparse


def replace_specific_variables(filename, targets):
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()

    out = []
    inside = None

    for line in lines:
        stripped = line.strip()
        if stripped in targets:
            inside = stripped
        if stripped == "}":
            inside = None
        if inside and "type" in stripped and "patch" in stripped:
            line = line.replace("patch", "wall")
        if inside and "physicalType" in stripped and "patch" in stripped:
            line = line.replace("patch", "wall")
        out.append(line)

    with open(filename, "w", encoding="utf-8") as f:
        f.writelines(out)


# CHANGE 'PATCH' TO 'WALL' FOR WALLS IN CONSTANT/POLYMESH/BOUNDARY

parser = argparse.ArgumentParser(description="Determines which case to run in.")

parser.add_argument("--breeder", type=str, help="Breeder material, 'lipb' or 'flibe'.")
parser.add_argument("--velocity", type=str, help="Velocity of OpenFOAM case being run.")

args = parser.parse_args()

case_velocity = args.velocity
breeder = args.breeder

replace_specific_variables(
    filename=breeder+"/vel_"
    + str(case_velocity)
    + "/constant/polyMesh/boundary",
    targets={"walls"},
)