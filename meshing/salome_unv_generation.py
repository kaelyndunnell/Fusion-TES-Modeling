import sys
import os
import salome
salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
import numpy as np
import GEOM
from salome.geom import geomBuilder
import SMESH
from salome.smesh import smeshBuilder


def build_and_export(tube_radius, bend_radius, length, unv_dir,
                     n_segments=100, growth_rate=1.01, radial_segments=8):
    """
    Build a U-tube hex mesh and export to UNV.

    Parameters
    ----------
    tube_radius      : float  - tube cross-section radius
    bend_radius      : float  - U-bend arc radius
    length           : float  - straight leg length
    unv_dir          : str    - output directory for .unv file
    n_segments       : int    - number of axial segments along path
    growth_rate      : float  - axial segment growth rate
    radial_segments  : int    - number of radial/circumferential segments (fixed)

    Returns
    -------
    str  - path to exported .unv file, or None on failure
    """
    if tube_radius >= bend_radius:
        print("[SKIP] tube_radius ({}) must be < bend_radius ({})".format(
            tube_radius, bend_radius))
        return None

    os.makedirs(unv_dir, exist_ok=True)

    geompy = geomBuilder.New()
    smesh  = smeshBuilder.New()

    # --- Geometry ---
    geomObj_1 = geompy.MakeMarker(0, 0, 0, 1, 0, 0, 0, 1, 0)
    sk = geompy.Sketcher2D()
    sk.addPoint(0.0, 0.0)
    sk.addSegmentLength(float(length))
    sk.addArcRadiusLength(float(bend_radius), 90.0)
    sk.addArcRadiusLength(float(bend_radius), 90.0)
    sk.addSegmentLength(float(length))
    Sketch_1 = sk.wire(geomObj_1)

    Divided_Disk_1 = geompy.MakeDividedDisk(float(tube_radius), 2, GEOM.SQUARE)
    Pipe_1 = geompy.MakePipe(Divided_Disk_1, Sketch_1)

    # --- Groups ---
    [Face_1,Face_2,Face_3,Face_4,Face_5,
     Face_6,Face_7,Face_8,Face_9,Face_10] = geompy.SubShapes(
        Pipe_1, [4, 220, 256, 94, 98, 157, 252, 161, 224, 258])
    [Solid_1,Solid_2,Solid_3,Solid_4,Solid_5] = geompy.ExtractShapes(
        Pipe_1, geompy.ShapeType["SOLID"], True)

    inlet  = geompy.CreateGroup(Pipe_1, geompy.ShapeType["FACE"])
    geompy.UnionIDs(inlet,  [4, 256, 98, 161, 224])
    outlet = geompy.CreateGroup(Pipe_1, geompy.ShapeType["FACE"])
    geompy.UnionIDs(outlet, [220, 94, 157, 252, 258])
    walls  = geompy.CreateGroup(Pipe_1, geompy.ShapeType["FACE"])
    geompy.UnionIDs(walls,  [52, 57, 105, 47, 120, 42, 115, 110,
                              173, 178, 234, 231, 168, 228, 237, 183])
    fluid  = geompy.CreateGroup(Pipe_1, geompy.ShapeType["SOLID"])
    geompy.UnionIDs(fluid,  [222, 2, 254, 159, 96])

    # --- Identify axial edges (swept path edges, much longer than radial edges) ---
    all_edges   = geompy.ExtractShapes(Pipe_1, geompy.ShapeType["EDGE"], True)

    axial_edges = []
    circumferential_edges = []

    for edge in all_edges:
        edge_length = geompy.BasicProperties(edge)[0]
        if edge_length > float(tube_radius) * 3:
            axial_edges.append(edge)
        else:
            if edge_length < float(tube_radius)/2: # edges from square edge to outer circumference of pipe
                circumferential_edges.append(edge)
    print("  Found {} axial edges".format(len(axial_edges)))
    print("  Found {} circumferential edges".format(len(circumferential_edges)))


    axial_group = geompy.CreateGroup(Pipe_1, geompy.ShapeType["EDGE"])
    geompy.UnionList(axial_group, axial_edges)

    circum_group = geompy.CreateGroup(Pipe_1, geompy.ShapeType["EDGE"])
    geompy.UnionList(circum_group, circumferential_edges)

    # --- Mesh ---
    mesh_name = "u_tube_{}_{}_{}" .format(tube_radius, bend_radius, length)
    Mesh_1 = smesh.Mesh(fluid, mesh_name)

    # Global hypothesis — controls radial/circumferential edges
    Regular_1D_global = Mesh_1.Segment()
    Regular_1D_global.NumberOfSegments(radial_segments, growth_rate, []) # number of segments along square edge in center
    Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
    Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)

    # Sub-mesh on axial edges — controls axial resolution independently
    Regular_1D_axial = Mesh_1.Segment(geom=axial_group)
    Regular_1D_axial.NumberOfSegments(n_segments, 1.0, []) # keep axial growth rate constant

    # sub-mesh on circumferential edges 
    Regular_1D_circum = Mesh_1.Segment(geom=circum_group)
    Regular_1D_circum.NumberOfSegments(20, 1.2, []) 

    isDone = Mesh_1.Compute()
    Mesh_1.CheckCompute()
    smesh.SetName(Mesh_1, mesh_name)

    smesh.SetName(Mesh_1.GroupOnGeom(inlet,  'inlet',  SMESH.FACE),   'inlet')
    smesh.SetName(Mesh_1.GroupOnGeom(outlet, 'outlet', SMESH.FACE),   'outlet')
    smesh.SetName(Mesh_1.GroupOnGeom(walls,  'walls',  SMESH.FACE),   'walls')
    smesh.SetName(Mesh_1.GroupOnGeom(fluid,  'fluid',  SMESH.VOLUME), 'fluid')

    # --- Export ---
    unv_path = os.path.join(unv_dir, "{}.unv".format(mesh_name))
    try:
        Mesh_1.ExportUNV(unv_path)
        print("[OK] UNV: {}".format(unv_path))
        return unv_path
    except Exception as e:
        print("[ERROR] UNV export failed: {}".format(e))
        return None


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    tube_r  = float(os.environ.get('TUBE_R',  0.0046))
    bend_r  = float(os.environ.get('BEND_R',  10))
    length  = float(os.environ.get('LENGTH',  40))
    unv_dir = os.environ.get('UNV_DIR')

    print(f"--- tube_radius={tube_r}, bend_radius={bend_r}, length={length} ---")
    try:
        build_and_export(tube_r, bend_r, length, unv_dir)
    except Exception as e:
        print("[ERROR] {},{},{}: {}".format(tube_r, bend_r, length, e))

    print("=== Done ===")

    if salome.sg.hasDesktop():
        salome.sg.updateObjBrowser()