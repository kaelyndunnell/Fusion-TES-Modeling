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
parser.add_argument(
    "--geometry", type=str, help="Geometry of the case, 'inlet_tank' or 'pipe'."
)
parser.add_argument(
    "--velocity", type=str, help="Inlet velocity of OpenFOAM case, in m/s."
)
parser.add_argument(
    "--bend_radius", type=str, help="Bend radius of OpenFOAM case, in m."
)
parser.add_argument(
    "--length", type=str, help="Length of OpenFOAM case, in m."
)

args = parser.parse_args()

geometry = args.geometry
breeder = args.breeder
velocity = args.velocity
bend_radius = args.bend_radius
length = args.length

replace_specific_variables(
    filename=breeder
    + "/"
    + geometry
    + "_"
    + velocity
    + "m_s_r"
    + bend_radius
    + "_l"
    + length
    + "/constant/polyMesh/boundary",
    targets={"walls", "tank_walls"},
)
