import gmsh
from collections import defaultdict

print(f"GMSH Version:", gmsh.__version__)

gmsh.initialize()
gmsh.option.setString("Geometry.OCCTargetUnit", "M")
gmsh.model.add("festim_mesh")

cad_file_path = "/Users/ckhurana/FESTIM/FESTIM-dev/openfoam/OpenFOAM/TES-PAV/meshing/shell_reverse_scaled1000_12_9.step"
entities = gmsh.model.occ.importShapes(cad_file_path)
gmsh.model.occ.synchronize()


# Extract raw volumes
volumes = gmsh.model.occ.getEntities(3)
# print(f"Extracted {len(volumes)} raw volumes from CAD.")

# Fragment volumes
# print("Fragmenting volumes...")
# gmsh.model.occ.fragment(volumes, [])
gmsh.model.occ.synchronize()

# Final volumes
final_volumes = gmsh.model.getEntities(3)
print(f"Final number of volumes: {len(final_volumes)}")
total_surfaces = [s[1] for s in gmsh.model.getEntities(2)]
print(f"Total number of surfaces: {len(total_surfaces)}")

###############################################
# Build volume→surface mapping
###############################################
volume_surfaces = {}

for dim, vol_tag in final_volumes:
    boundary = gmsh.model.getBoundary([(3, vol_tag)], oriented=False, combined=False)
    surf_tags = [s[1] for s in boundary if s[0] == 2]
    volume_surfaces[vol_tag] = surf_tags

# Count occurrences to find shared surfaces
surface_counts = defaultdict(int)
for surf_list in volume_surfaces.values():
    for s in surf_list:
        surface_counts[s] += 1

shared_surfaces = [s for s, c in surface_counts.items() if c > 1]
print("Shared surfaces:", shared_surfaces)

###############################################
# Identify fluid vs pipes
###############################################
fluid_volume = 1  
pipe_volumes = [v[1] for v in final_volumes if v[1] != fluid_volume]

###############################################
# INLET / OUTLET (via curve loops)
###############################################
inlet_id = [42]
outlet_id = [41]

gmsh.model.occ.synchronize()

###############################################
# VACUUM SURFACES
###############################################
vacuum_surfaces = []

for v in pipe_volumes:
    for s in volume_surfaces[v]:
        if (s not in shared_surfaces):
            vacuum_surfaces.append(s)

print("Vacuum surfaces:", vacuum_surfaces)
walls_surfaces = [
    s for s in total_surfaces
    if s not in vacuum_surfaces
    and s not in (inlet_id, outlet_id)
]
print("Walls:", walls_surfaces)

###############################################
# PHYSICAL GROUPS
###############################################
gmsh.model.addPhysicalGroup(3, [fluid_volume], name="fluid")
gmsh.model.addPhysicalGroup(3, pipe_volumes, name="pipes")

gmsh.model.addPhysicalGroup(2, inlet_id, name="inlet")
gmsh.model.addPhysicalGroup(2, outlet_id, name="outlet")

gmsh.model.addPhysicalGroup(2, walls_surfaces, name="walls")
gmsh.model.addPhysicalGroup(2, vacuum_surfaces, name="vacuum")

###############################################
# MESH REFINEMENT
###############################################

# Global settings
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
# gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.001)
# gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 50)

# Define mesh sizes
pipe_mesh_size = 10.0      # Coarse mesh in pipes
fluid_mesh_size = 15      # Fine mesh in fluid
wall_mesh_size = 1      # Refined at walls
inlet_outlet_size = 1.0    # Very fine at inlet/outlet

# Set mesh size for PIPE volumes (coarse)
for vol in pipe_volumes:
    # Get all points in this volume
    volume_boundary = gmsh.model.getBoundary([(3, vol)], recursive=True)
    points = [p[1] for p in volume_boundary if p[0] == 0]
    for pt in points:
        gmsh.model.mesh.setSize([(0, pt)], pipe_mesh_size)

# Set mesh size for FLUID volume (fine)
fluid_boundary = gmsh.model.getBoundary([(3, fluid_volume)], recursive=True)
fluid_points = [p[1] for p in fluid_boundary if p[0] == 0]
for pt in fluid_points:
    gmsh.model.mesh.setSize([(0, pt)], fluid_mesh_size)

# WALL REFINEMENT (shared surfaces between fluid and pipes)
for s in shared_surfaces:
    # Get curves on this surface
    curves = gmsh.model.getBoundary([(2, s)], oriented=False)
    for (dim, cid) in curves:
        if dim == 1:
            # Get points on this curve
            curve_points = gmsh.model.getBoundary([(1, cid)], oriented=False)
            for (pdim, pid) in curve_points:
                if pdim == 0:
                    gmsh.model.mesh.setSize([(0, pid)], wall_mesh_size)

# INLET REFINEMENT
inlet_tag = inlet_id[0]
inlet_curves = gmsh.model.getBoundary([(2, inlet_tag)], oriented=False)
for (dim, cid) in inlet_curves:
    if dim == 1:
        inlet_points = gmsh.model.getBoundary([(1, cid)], oriented=False)
        for (pdim, pid) in inlet_points:
            if pdim == 0:
                gmsh.model.mesh.setSize([(0, pid)], inlet_outlet_size)

# OUTLET REFINEMENT
outlet_tag = outlet_id[0]
outlet_curves = gmsh.model.getBoundary([(2, outlet_tag)], oriented=False)
for (dim, cid) in outlet_curves:
    if dim == 1:
        outlet_points = gmsh.model.getBoundary([(1, cid)], oriented=False)
        for (pdim, pid) in outlet_points:
            if pdim == 0:
                gmsh.model.mesh.setSize([(0, pid)], inlet_outlet_size)

# Optional: Use mesh size fields for more control
# This creates a distance-based refinement from walls
use_distance_field = True
if use_distance_field:
    # Create distance field from walls
    field_id = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(field_id, "SurfacesList", shared_surfaces)
    
    # Create threshold field (fine near walls, coarse away)
    threshold_id = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold_id, "InField", field_id)
    gmsh.model.mesh.field.setNumber(threshold_id, "SizeMin", wall_mesh_size)
    gmsh.model.mesh.field.setNumber(threshold_id, "SizeMax", fluid_mesh_size)
    gmsh.model.mesh.field.setNumber(threshold_id, "DistMin", 0.5)  # Distance from wall where refinement starts
    gmsh.model.mesh.field.setNumber(threshold_id, "DistMax", 5.0)  # Distance where coarse mesh begins
    
    # Set as background field
    gmsh.model.mesh.field.setAsBackgroundMesh(threshold_id)

###############################################
# MESH GENERATION OPTIONS
# ###############################################
# gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay for 2D
# gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT for 3D (fast, good quality)
gmsh.option.setNumber("Mesh.Optimize", 1)  # Optimize mesh quality
gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)  # Additional optimization

###############################################
# GENERATE MESH
###############################################
# gmsh.option.
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.model.mesh.generate(3)
gmsh.write("coarser_multi_region_12_17.msh")
gmsh.fltk.run()
gmsh.finalize()
