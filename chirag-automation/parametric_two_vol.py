import subprocess
import os
import itertools

SCRIPT_DIR      = '/home/ckhurana/Fusion-TES-Modeling/chirag-automation'
SALOME_PATH     = '/opt/salome/salome'
MESH_SCRIPT     = os.path.join(SCRIPT_DIR, 'two_vol_salome.py')
UNV_DIR         = os.path.join(SCRIPT_DIR, 'meshes/bl/two_vol/UNV')

# ---------------------------------------------------------------------------
# Parameter sweep — edit these lists
# ---------------------------------------------------------------------------
TUBE_RADII   = [0.00465]          # m
BL_THICK     = [r / 8 for r in TUBE_RADII]   # m
WALL_THICK   = [0.0004]         # m
BEND_RADII   = [0.5, 2.5]  # m
LENGTHS      = [1.0, 8.0]   # m

N_SEGMENTS   = 200
GROWTH_RATE  = 1.0
RADIAL_SEGS  = 8
# ---------------------------------------------------------------------------


def mesh_utube(tube_r, bl_t, wall_t, bend_r, length, unv_dir):
    env = os.environ.copy()
    env['TUBE_R']      = str(tube_r)
    env['BL_T']        = str(bl_t)
    env['WALL_T']      = str(wall_t)
    env['BEND_R']      = str(bend_r)
    env['LENGTH']      = str(length)
    env['UNV_DIR']     = unv_dir
    env['N_SEGMENTS']  = str(N_SEGMENTS)
    env['GROWTH_RATE'] = str(GROWTH_RATE)
    env['RADIAL_SEGS'] = str(RADIAL_SEGS)
    env.pop('LD_LIBRARY_PATH', None)
    env.pop('CONDA_PREFIX', None)
    env.pop('CONDA_DEFAULT_ENV', None)
    env['PYTHONWARNINGS'] = 'ignore::SyntaxWarning'

    print(f"Meshing tube_r={tube_r} bl_t={bl_t} wall_t={wall_t} "
          f"bend_r={bend_r} length={length}...")

    result = subprocess.run(
        [SALOME_PATH, '-t', 'python', MESH_SCRIPT],
        env=env, capture_output=True, text=True,
        cwd='/opt/salome'
    )
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(f"Meshing failed: {tube_r},{bl_t},{wall_t},{bend_r},{length}")

    mesh_path = os.path.join(
        unv_dir,
        f"u_tube_{float(tube_r)}_{float(bl_t)}_{float(wall_t)}_"
        f"{float(bend_r)}_{float(length)}.unv"
    )
    if not os.path.exists(mesh_path):
        raise RuntimeError(f"UNV not created: {mesh_path}")
    return mesh_path


def run_fem(mesh_path, **params):
    print(f"Running FEM on {mesh_path}...")
    pass


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------
os.makedirs(UNV_DIR, exist_ok=True)

params = list(itertools.product(
    TUBE_RADII, BL_THICK, WALL_THICK, BEND_RADII, LENGTHS
))

print(f"Total cases: {len(params)}")

for tube_r, bl_t, wall_t, bend_r, length in params:
    # skip invalid geometry
    if tube_r >= bend_r:
        print(f"[SKIP] tube_r={tube_r} >= bend_r={bend_r}")
        continue
    if tube_r - bl_t <= 0:
        print(f"[SKIP] bl_t={bl_t} >= tube_r={tube_r}")
        continue
    try:
        mesh_file = mesh_utube(tube_r, bl_t, wall_t, bend_r, length, UNV_DIR)
        run_fem(mesh_file, tube_r=tube_r, bl_t=bl_t, wall_t=wall_t,
                bend_r=bend_r, length=length)
    except RuntimeError as e:
        print(f"[ERROR] {e}")

print("\n=== All done ===")