import sys
import os
import math
import salome
salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()

import GEOM
from salome.geom import geomBuilder
import SMESH
from salome.smesh import smeshBuilder


def build_and_export(tube_radius, bl_thickness, wall_thickness, bend_radius, length, med_dir,
                     n_segments=60, growth_rate=1.0, radial_segments=15):
    """
    Build a U-tube hex mesh with fluid + BL + wall regions and export to UNV.
    Uses hardcoded subshape IDs verified in Salome GUI.
    """
    r1 = float(tube_radius) - float(bl_thickness)
    r2 = float(tube_radius)
    r3 = r2 + float(wall_thickness)

    if r1 >= float(bend_radius):
        print("[SKIP] tube_radius ({}) must be < bend_radius ({})".format(r1, bend_radius))
        return None

    # os.makedirs(unv_dir, exist_ok=True)
    os.makedirs(med_dir, exist_ok=True)  # add this

    geompy = geomBuilder.New()
    smesh  = smeshBuilder.New()

    # --- Geometry: sketch path ---
    geomObj_1 = geompy.MakeMarker(0, 0, 0, 1, 0, 0, 0, 1, 0)
    sk = geompy.Sketcher2D()
    sk.addPoint(0.0, 0.0)
    sk.addSegmentLength(float(length))
    sk.addArcRadiusLength(float(bend_radius), 90.0)
    sk.addArcRadiusLength(float(bend_radius), 90.0)
    sk.addSegmentLength(float(length))
    Sketch_1 = sk.wire(geomObj_1)

    # --- Three concentric divided disks ---
    Divided_Disk_1 = geompy.MakeDividedDisk(r1, 2, GEOM.SQUARE)
    Divided_Disk_2 = geompy.MakeDividedDisk(r2, 2, GEOM.SQUARE)
    Divided_Disk_3 = geompy.MakeDividedDisk(r3, 2, GEOM.SQUARE)

    # --- Sweep along path ---
    Pipe_1 = geompy.MakePipe(Divided_Disk_1, Sketch_1)
    Pipe_2 = geompy.MakePipe(Divided_Disk_2, Sketch_1)
    Pipe_3 = geompy.MakePipe(Divided_Disk_3, Sketch_1)

    # --- Two cuts ---
    Cut_1 = geompy.MakeCutList(Pipe_3, [Pipe_2], True)
    Cut_2 = geompy.MakeCutList(Pipe_2, [Pipe_1], True)

    # --- Partition ---
    two_vol = geompy.MakePartition(
        [Pipe_1, Cut_1, Cut_2], [], [], [],
        geompy.ShapeType["SOLID"], 0, [], 0
    )

    # --- Groups using hardcoded IDs verified in Salome GUI ---
    fluid = geompy.CreateGroup(two_vol, geompy.ShapeType["SOLID"])
    geompy.UnionIDs(fluid, [582, 512, 222, 2, 254, 159, 96, 563, 544])

    wall = geompy.CreateGroup(two_vol, geompy.ShapeType["SOLID"])
    geompy.UnionIDs(wall, [480, 260, 417, 354])

    inlet = geompy.CreateGroup(two_vol, geompy.ShapeType["FACE"])
    geompy.UnionIDs(inlet, [252, 157, 561, 94, 258, 586, 220, 580, 542])

    outlet = geompy.CreateGroup(two_vol, geompy.ShapeType["FACE"])
    geompy.UnionIDs(outlet, [565, 546, 161, 98, 256, 224, 4, 584, 514])

    vacuum = geompy.CreateGroup(two_vol, geompy.ShapeType["FACE"])
    geompy.UnionIDs(vacuum, [373, 489, 284, 436, 507,
                              344, 324, 386, 449, 399, 462, 475, 495, 304,
                              412, 501])

    walls = geompy.CreateGroup(two_vol, geompy.ShapeType["FACE"])
    geompy.UnionIDs(walls, [510, 415, 352, 478, 419, 482, 262, 356])

    interfaces = geompy.CreateGroup(two_vol, geompy.ShapeType["FACE"])
    geompy.UnionIDs(interfaces, [504, 465, 402, 332, 452, 389, 498, 312, 492, 292, 439, 376, 486, 272, 426, 363])

    # --- Identify axial edges for sub-mesh ---
    all_edges   = geompy.ExtractShapes(two_vol, geompy.ShapeType["EDGE"], True)
    axial_edges = [e for e in all_edges
                   if geompy.BasicProperties(e)[0] > r3 * 2]
    print("  Axial edges: {}".format(len(axial_edges)))

    axial_group = geompy.CreateGroup(two_vol, geompy.ShapeType["EDGE"])
    geompy.UnionList(axial_group, axial_edges)

    # --- Mesh ---
    mesh_name = "two_vol_{}_{:.2f}_{:.2f}".format(
        tube_radius, bend_radius, length)
    Mesh_1 = smesh.Mesh(two_vol, mesh_name)

    # Global hypothesis
    Regular_1D_global = Mesh_1.Segment()
    Regular_1D_global.NumberOfSegments(radial_segments, 1.0, [])
    Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
    Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)

    # Axial sub-mesh
    Regular_1D_axial = Mesh_1.Segment(geom=axial_group)
    Regular_1D_axial.NumberOfSegments(n_segments, growth_rate, [])

    isDone = Mesh_1.Compute()
    Mesh_1.CheckCompute()
    smesh.SetName(Mesh_1, mesh_name)

    smesh.SetName(Mesh_1.GroupOnGeom(fluid,  'fluid',  SMESH.VOLUME), 'fluid')
    smesh.SetName(Mesh_1.GroupOnGeom(wall,   'membrane',   SMESH.VOLUME), 'membrane')
    smesh.SetName(Mesh_1.GroupOnGeom(inlet,  'inlet',  SMESH.FACE),   'inlet')
    smesh.SetName(Mesh_1.GroupOnGeom(outlet, 'outlet', SMESH.FACE),   'outlet')
    smesh.SetName(Mesh_1.GroupOnGeom(walls,  'walls',  SMESH.FACE),   'walls')
    smesh.SetName(Mesh_1.GroupOnGeom(vacuum,  'vacuum',  SMESH.FACE),   'vacuum')
    smesh.SetName(Mesh_1.GroupOnGeom(interfaces,  'interfaces',  SMESH.FACE),   'interfaces')

    # --- Export ---
    # unv_path = os.path.join(unv_dir, "{}.unv".format(mesh_name))
    med_path = os.path.join(med_dir, "{}.med".format(mesh_name))
    try:
        # Mesh_1.ExportUNV(unv_path)
        Mesh_1.ExportMED(med_path, 0, 41, 1, Mesh_1, 1, [], '', -1, 1)
        print("[OK] MED: {}".format(med_path))
        return med_path #, unv_path
    except Exception as e:
        print("[ERROR] UNV export failed: {}".format(e))
        return None, None


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    tube_r  = float(os.environ.get('TUBE_R',    4.6e-3))
    bl_t    = float(os.environ.get('BL_T',      tube_r/4))
    wall_t  = float(os.environ.get('WALL_T',    0.4e-3))
    bend_r  = float(os.environ.get('BEND_R',   0.5))
    length  = float(os.environ.get('LENGTH',   1.0))
    # unv_dir = os.environ.get('UNV_DIR')
    med_dir = os.environ.get('MESH_DIR')
    n_seg   = int(float(os.environ.get('N_SEGMENTS', round(length*65))))
    growth  = float(os.environ.get('GROWTH_RATE', 1.1))
    radial  = int(float(os.environ.get('RADIAL_SEGS', 8)))

    print("--- tube_r={} bl_t={} wall_t={} bend_r={} length={} ---".format(
        tube_r, bl_t, wall_t, bend_r, length))
    try:
        build_and_export(tube_r, bl_t, wall_t, bend_r, length, med_dir,
                         n_segments=n_seg,
                         growth_rate=growth,
                         radial_segments=radial)
    except Exception as e:
        print("[ERROR] {},{},{},{},{}: {}".format(
            tube_r, bl_t, wall_t, bend_r, length, e))

    print("=== Done ===")

    if salome.sg.hasDesktop():
        salome.sg.updateObjBrowser()