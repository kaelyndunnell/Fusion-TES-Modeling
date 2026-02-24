import gmsh

###############################################
###### CREATE FESTIM MESH FROM CAD MODEL ######
###############################################

# LOAD CAD AND INITIALIZE MESH
gmsh.initialize()
gmsh.option.setString(
    "Geometry.OCCTargetUnit", "M"
)  # make sure gmsh reads .step file in m
gmsh.model.add("tes-model-festim")  # reading cadquery file as if in M

cad_file_path = "meshing/CAD/tes_festim.brep"
gmsh.merge(cad_file_path)
gmsh.model.occ.synchronize()

# EXTRACT ALL VOLUMES
volumes = [e for e in gmsh.model.occ.getEntities() if e[0] == 3]
print(f"Extracted {len(volumes)} raw volumes from CAD.")
gmsh.model.occ.synchronize()

surfaces = gmsh.model.occ.getEntities(dim=2)

##### MESH SETTINGS #####

# Mesh settings
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.Algorithm3D", 4)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 15)

# Set reasonable mesh size limits to avoid issues with small features
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 1e-5)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.5)

# CRITICAL: Turn off compound entities and periodic matching
gmsh.option.setNumber("Mesh.CompoundClassify", 0)
gmsh.option.setNumber("Mesh.CompoundMeshSizeFactor", 1)
gmsh.model.occ.synchronize()

##### TAG & NAME PHYSICAL GROUPS #####
# need to open mesh in gmsh gui to determine the proper tagging as below

# volumes
breeder_marker = 1
permeable_membrane_marker = 2

gmsh.model.addPhysicalGroup(3, [1], breeder_marker, name=f"breeder")
gmsh.model.addPhysicalGroup(
    3,
    [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
    permeable_membrane_marker,
    name=f"permeable_membrane",
)

# interface surfaces
surfaces = gmsh.model.getEntities(dim=2)
interface_surfaces = [
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
interface_marker = 3
gmsh.model.addPhysicalGroup(2, interface_surfaces, interface_marker, name="interfaces")

# other surfaces
inlet = [34]
inlet_marker = 4

outlet = [72]
outlet_marker = 5

vacuum_surfaces = list(range(73, 163))
vacuum_marker = 6

wall_surfaces = [1, 2, 3, 4, 14, 24, 34, 35, 36, 37, 47, 66, 67, 68, 69, 70, 71, 72]
walls_marker = 7

gmsh.model.addPhysicalGroup(2, inlet, inlet_marker, name="inlet")
gmsh.model.addPhysicalGroup(2, outlet, outlet_marker, name="outlet")
gmsh.model.addPhysicalGroup(2, vacuum_surfaces, vacuum_marker, name="vacuum")
gmsh.model.addPhysicalGroup(2, wall_surfaces, walls_marker, name="walls")


##### MESH SIZE & REFINEMENT #####

gmsh.model.occ.synchronize()

# extra refinement for tanks
distance_field = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(distance_field, "FacesList", wall_surfaces)

# use threshold field to refine near surfaces
threshold_field = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold_field, "IField", distance_field)
gmsh.model.mesh.field.setNumber(
    threshold_field, "SizeMin", 1e-1
)  # smallest mesh size near surfaces
gmsh.model.mesh.field.setNumber(
    threshold_field, "SizeMax", 5e-1
)  # mesh size far from surfaces
gmsh.model.mesh.field.setNumber(
    threshold_field, "DistMin", 1.5e-1
)  # distance where within which mesh is fully refined
gmsh.model.mesh.field.setNumber(
    threshold_field, "DistMax", 2e-1
)  # distance where mesh transitions to coarse size

# set threshold field as background field (doesn't impact probe surface)
gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)

##### SYNC & GENERATE MESH #####
gmsh.model.occ.synchronize()
# gmsh.fltk.run()
gmsh.model.mesh.generate(3)

# # ##### SAVE MESH #####co
output_file = "meshing/tes_festim.msh"
gmsh.write(output_file)
gmsh.finalize()
