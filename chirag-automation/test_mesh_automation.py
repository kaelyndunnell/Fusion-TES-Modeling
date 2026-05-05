import subprocess
import os

SALOME_PATH      = '/opt/salome/salome'
MESH_SCRIPT_PATH = '/home/ckhurana/Fusion-TES-Modeling/chirag-automation/salome_unv_generation.py'


def mesh_utube(tube_r, bend_r, length, unv_dir):
    env = os.environ.copy()
    env['TUBE_R']  = str(tube_r)
    env['BEND_R']  = str(bend_r)
    env['LENGTH']  = str(length)
    env['UNV_DIR'] = unv_dir

    print(f"Meshing u_tube_{tube_r}_{bend_r}_{length}...")
    result = subprocess.run(
        [SALOME_PATH, '-t', 'python', MESH_SCRIPT_PATH],
        env=env, capture_output=True, text=True,
        cwd='/opt/salome'  # run from salome's own directory
    )
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(f"Meshing failed for {tube_r},{bend_r},{length}")

    # match the float filename that Salome actually produces
    mesh_path = os.path.join(unv_dir, f"u_tube_{float(tube_r)}_{float(bend_r)}_{float(length)}.unv")
    if not os.path.exists(mesh_path):
        raise RuntimeError(f"UNV file not created: {mesh_path}")
    return mesh_path


def run_fem(mesh_path, tube_r, bend_r, length):
    print(f"Running FEM on {mesh_path}...")
    pass


# --- Main sweep ---

# Define directory to store UNV files
UNV_DIR = '/home/ckhurana/Fusion-TES-Modeling/chirag-automation/meshes/UNV'

# Define parameter ranges (as lists)

for tube_r in [2]:
    for bend_r in [10, 20]:
        for length in [40, 80]:
            try:
                mesh_file = mesh_utube(tube_r, bend_r, length, UNV_DIR)
                run_fem(mesh_file, tube_r, bend_r, length)
            except RuntimeError as e:
                print(f"[ERROR] {e}")
