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

print("Fragmenting volumes to define interfaces...")
gmsh.model.occ.fragment(volumes, [])
gmsh.model.occ.synchronize()

# FINAL VOLUMES AFTER FRAGMENT
final_volumes = gmsh.model.getEntities(dim=3)
print(f"Final number of volumes: {len(final_volumes)}")

##### MESH SETTINGS #####

# Mesh settings
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)

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
interface_surfaces = list(range(5, 74)) + [
    76,
    78,
    81,
    83,
    86,
    88,
    91,
    93,
    96,
    98,
    101,
    103,
    106,
    108,
    111,
    113,
    116,
    118,
    121,
    123,
    126,
    128,
    131,
    133,
    136,
    138,
    141,
    143,
    146,
    148,
    151,
    153,
    156,
    158,
    161,
]
interface_marker = 3
gmsh.model.addPhysicalGroup(2, interface_surfaces, interface_marker, name="interfaces")

# other surfaces
inlet = [169]
inlet_marker = 4

outlet = [180]
outlet_marker = 5

vacuum_surfaces = [
    74,
    75,
    77,
    79,
    80,
    82,
    84,
    85,
    87,
    89,
    90,
    92,
    94,
    95,
    97,
    99,
    100,
    102,
    104,
    105,
    107,
    109,
    110,
    112,
    114,
    115,
    117,
    119,
    120,
    122,
    124,
    125,
    127,
    129,
    130,
    132,
    134,
    135,
    137,
    139,
    140,
    142,
    144,
    145,
    147,
    149,
    150,
    152,
    154,
    155,
    157,
    159,
    160,
]
vacuum_marker = 6

wall_surfaces = [
    163,
    164,
    165,
    166,
    167,
    168,
    170,
    171,
    172,
    173,
    174,
    175,
    176,
    177,
    178,
    179,
]
walls_marker = 7

gmsh.model.addPhysicalGroup(2, inlet, inlet_marker, name="inlet")
gmsh.model.addPhysicalGroup(2, outlet, outlet_marker, name="outlet")
gmsh.model.addPhysicalGroup(2, vacuum_surfaces, vacuum_marker, name="vacuum")
gmsh.model.addPhysicalGroup(2, wall_surfaces, walls_marker, name="walls")


##### MESH SIZE & REFINEMENT #####

gmsh.model.occ.synchronize()

# # extra refinement for tanks --> TODO
# distance_field = gmsh.model.mesh.field.add("Distance")
# gmsh.model.mesh.field.setNumbers(distance_field, "FacesList", wall_surfaces)

# # use threshold field to refine near surfaces
# threshold_field = gmsh.model.mesh.field.add("Threshold")
# gmsh.model.mesh.field.setNumber(threshold_field, "IField", distance_field)
# gmsh.model.mesh.field.setNumber(
#     threshold_field, "SizeMin", 1e-1
# )  # smallest mesh size near surfaces
# gmsh.model.mesh.field.setNumber(
#     threshold_field, "SizeMax", 5e-1
# )  # mesh size far from surfaces
# gmsh.model.mesh.field.setNumber(
#     threshold_field, "DistMin", 1.5e-1
# )  # distance where within which mesh is fully refined
# gmsh.model.mesh.field.setNumber(
#     threshold_field, "DistMax", 2e-1
# )  # distance where mesh transitions to coarse size

# # set threshold field as background field (doesn't impact probe surface)
# gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)

##### SYNC & GENERATE MESH #####
gmsh.model.occ.synchronize()
# gmsh.fltk.run()
gmsh.model.mesh.generate(3)

# # ##### SAVE MESH #####co
output_file = "meshing/tes_festim.msh"
gmsh.write(output_file)
gmsh.finalize()
