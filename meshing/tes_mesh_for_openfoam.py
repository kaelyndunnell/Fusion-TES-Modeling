import gmsh

##### LOAD CAD AND INITIALIZE MESH #####

gmsh.initialize()
gmsh.option.setString(
    "Geometry.OCCTargetUnit", "MM"
)  # make sure gmsh reads .step file in mm, what CADQuery exports in
gmsh.model.add("tes-model")

cad_file_path = "meshing/CAD/tes_openfoam.step"

entities = gmsh.model.occ.importShapes(cad_file_path)
gmsh.model.occ.synchronize()

##### TAG & NAME PHYSICAL GROUPS #####
# need to open mesh in gmsh gui to determine the proper tagging as below

volumes = [e for e in gmsh.model.occ.getEntities() if e[0] == 3]
print(f"Extracted {len(volumes)} raw volumes from CAD.")

# volumes
vol_marker = 1  # fluid volume

gmsh.model.addPhysicalGroup(3, [1], vol_marker, name="fluid")

# surfaces
surfaces = gmsh.model.getEntities(dim=2)

inlet_tank = [2, 3, 4, 14, 24]
inlet = [34]
outlet_tank = [67, 68, 69, 70, 71]
outlet = [72]
mixing_tank = [
    36,
    37,
    47,
]
vacuum = [
    1,
    66,
    35,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    23,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    32,
    33,
    38,
    39,
    40,
    41,
    42,
    43,
    44,
    45,
    46,
    48,
    49,
    50,
    51,
    52,
    53,
    54,
    55,
    56,
    57,
    58,
    59,
    60,
    61,
    62,
    63,
    64,
    65,
]

inlet_tank_marker = 2
inlet_marker = 3
outlet_tank_marker = 4
outlet_marker = 5
mixing_tank_marker = 6
vacuum_marker = 7

gmsh.model.addPhysicalGroup(2, inlet_tank, inlet_tank_marker, name="inlet_tank")
gmsh.model.addPhysicalGroup(2, inlet, inlet_marker, name="inlet")
gmsh.model.addPhysicalGroup(2, outlet_tank, outlet_tank_marker, name="outlet_tank")
gmsh.model.addPhysicalGroup(2, outlet, outlet_marker, name="outlet")
gmsh.model.addPhysicalGroup(2, mixing_tank, mixing_tank_marker, name="mixing_tank")
gmsh.model.addPhysicalGroup(2, vacuum, vacuum_marker, name="vacuum")

gmsh.model.occ.synchronize()

##### MESH SIZE & REFINEMENT #####

# want finer mesh size at the walls
walls = (
    inlet_tank + outlet_tank + mixing_tank + vacuum
)  # use the internal gmsh labels, NOT the markers!

distance_field = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(distance_field, "FacesList", walls)

threshold_field = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold_field, "IField", distance_field)
gmsh.model.mesh.field.setNumber(
    threshold_field, "SizeMin", 0.05
)  # smallest mesh size near surfaces
gmsh.model.mesh.field.setNumber(
    threshold_field, "SizeMax", 0.1
)  # mesh size far from surfaces
gmsh.model.mesh.field.setNumber(
    threshold_field, "DistMin", 0.15
)  # distance where within which mesh is fully refined
gmsh.model.mesh.field.setNumber(
    threshold_field, "DistMax", 0.2
)  # distance where mesh transitions to coarse size

# set threshold field as background field
gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)

##### SYNC & GENERATE MESH #####
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)

##### SAVE MESH #####
output_file = "meshing/tes_openfoam.msh"
gmsh.write(output_file)
gmsh.finalize()
