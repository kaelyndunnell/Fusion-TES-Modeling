
###
### Parametric U-tube hex mesh generator for SALOME v9.14.0
###
### Three parameters: tube radius (disk), bend radius (arc), leg length
###
### Exports to:
###   parametric/u_tube/MED/u_tube_[tube_r]_[bend_r]_[length].med    (solver use)
###   parametric/u_tube/XDMF/u_tube_[tube_r]_[bend_r]_[length].xdmf  (ParaView)
###
### Groups inlet / outlet / walls / fluid are visible in both formats.
###
### First run only: SALOME needs meshio installed.
### In the SALOME terminal run:  pip install meshio --break-system-packages
###

import sys
import os
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()

# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------

TUBE_RADII  = [2, 4]        # tube cross-section radius (disk radius)
BEND_RADII  = [10, 20]      # U-bend arc radius in sketch
LENGTHS     = [40, 80]      # straight-leg length on each side of the U

# NOTE: TUBE_RADIUS must be < BEND_RADIUS or geometry will self-intersect

N_SEGMENTS  = 30
GROWTH_RATE = 1.3

# ---------------------------------------------------------------------------
# Output folders
# ---------------------------------------------------------------------------
SCRIPT_DIR = r'C:/Users/ckhur/Desktop/PSFC/FESTIM-DEV/TES/SALOME'
BASE_OUT   = os.path.join(SCRIPT_DIR, "testing_groups")
MED_DIR    = os.path.join(BASE_OUT, "MED")
XDMF_DIR   = os.path.join(BASE_OUT, "XDMF")

os.makedirs(MED_DIR,  exist_ok=True)
os.makedirs(XDMF_DIR, exist_ok=True)

# ---------------------------------------------------------------------------

import GEOM
from salome.geom import geomBuilder
import SMESH
from salome.smesh import smeshBuilder


def build_and_export(tube_radius, bend_radius, length):
    """Build one U-tube, assign boundary groups, export MED + XDMF."""

    if tube_radius >= bend_radius:
        print("[SKIP] tube_radius ({}) must be < bend_radius ({})".format(
            tube_radius, bend_radius))
        return

    geompy = geomBuilder.New()
    smesh  = smeshBuilder.New()

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # GEOM boundary groups
    # ------------------------------------------------------------------
    # all_faces  = geompy.ExtractShapes(Pipe_1, geompy.ShapeType["FACE"],  True)
    # all_solids = geompy.ExtractShapes(Pipe_1, geompy.ShapeType["SOLID"], True)

    # sk_verts  = geompy.ExtractShapes(Sketch_1, geompy.ShapeType["VERTEX"], True)
    # pt_start  = sk_verts[0]
    # pt_end    = sk_verts[-1]
    # threshold = float(tube_radius) * 1.5

    [Face_1,Face_2,Face_3,Face_4,Face_5,Face_6,Face_7,Face_8,Face_9,Face_10] = geompy.SubShapes(Pipe_1, [4, 220, 256, 94, 98, 157, 252, 161, 224, 258])
    [Solid_1,Solid_2,Solid_3,Solid_4,Solid_5] = geompy.ExtractShapes(Pipe_1, geompy.ShapeType["SOLID"], True)
    inlet = geompy.CreateGroup(Pipe_1, geompy.ShapeType["FACE"])
    geompy.UnionIDs(inlet, [4, 256, 98, 161, 224])
    outlet = geompy.CreateGroup(Pipe_1, geompy.ShapeType["FACE"])
    geompy.UnionIDs(outlet, [220, 94, 157, 252, 258])
    walls = geompy.CreateGroup(Pipe_1, geompy.ShapeType["FACE"])
    geompy.UnionIDs(walls, [52, 57, 105, 47, 120, 42, 115, 110, 173, 178, 234, 231, 168, 228, 237, 183])
    fluid = geompy.CreateGroup(Pipe_1, geompy.ShapeType["SOLID"])
    geompy.UnionIDs(fluid, [222, 2, 254, 159, 96])
    
    # sk_verts = geompy.ExtractShapes(Sketch_1, geompy.ShapeType["VERTEX"], True)
    # print("Sketch vertices:", len(sk_verts))
    # for i, v in enumerate(sk_verts):
    #     print(f"  {i}: {geompy.PointCoordinates(v)}")

    # print(f"\ntube_radius={tube_radius}, threshold={tube_radius * 1.5}")
    # print("Face classifications:")
    # for face in all_faces:
    #     fid = geompy.GetSubShapeID(Pipe_1, face)
    #     d0  = geompy.MinDistance(face, sk_verts[0])
    #     d1  = geompy.MinDistance(face, sk_verts[-1])
    #     area = geompy.BasicProperties(face)[1]
    #     print(f"  face_id={fid:4d}  d_start={d0:8.3f}  d_end={d1:8.3f}  area={area:.3f}")
    # ------------------------------------------------------------------
    # Meshing
    # ------------------------------------------------------------------
    mesh_name = "u_tube_{}_{}_{}" .format(tube_radius, bend_radius, length)
    Mesh_1 = smesh.Mesh(fluid, mesh_name)

    Regular_1D = Mesh_1.Segment()
    Regular_1D.NumberOfSegments(N_SEGMENTS, GROWTH_RATE, [])
    Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
    Mesh_1.Hexahedron(algo=smeshBuilder.Hexa)

    isDone = Mesh_1.Compute()
    Mesh_1.CheckCompute()
    smesh.SetName(Mesh_1, mesh_name)

    # ------------------------------------------------------------------
    # Transfer GEOM groups -> SMESH groups (embedded in exported files)
    # ------------------------------------------------------------------
    smesh.SetName(Mesh_1.GroupOnGeom(inlet,  'inlet',  SMESH.FACE),   'inlet')
    smesh.SetName(Mesh_1.GroupOnGeom(outlet, 'outlet', SMESH.FACE),   'outlet')
    smesh.SetName(Mesh_1.GroupOnGeom(walls,  'walls',  SMESH.FACE),   'walls')
    smesh.SetName(Mesh_1.GroupOnGeom(fluid,  'fluid',  SMESH.VOLUME), 'fluid')

    # ------------------------------------------------------------------
    # Export MED
    # ------------------------------------------------------------------
    med_path = os.path.join(MED_DIR, "{}.med".format(mesh_name))
    try:
        Mesh_1.ExportMED(med_path, 1, 41, 1, Mesh_1, 1, [], '', -1, 1)
        print("[OK] MED: {}".format(med_path))
    except Exception as e:
        print("[ERROR] MED export failed: {}".format(e))
        return

# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------
for tube_r in TUBE_RADII:
    for bend_r in BEND_RADII:
        for length in LENGTHS:
            print("\n--- tube_radius={}, bend_radius={}, length={} ---".format(
                tube_r, bend_r, length))
            try:
                build_and_export(tube_r, bend_r, length)
            except Exception as e:
                print("[ERROR] {},{},{}: {}".format(tube_r, bend_r, length, e))

print("\n=== Done ===")
print("  MED  (solvers)  -> {}".format(MED_DIR))

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()