#!/usr/bin/env python
import os, sys
import numpy as np
import meshio

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MED_DIR    = os.path.join(SCRIPT_DIR, "meshes", "MED")
MSH_DIR    = os.path.join(SCRIPT_DIR, "meshes", "MSH")
os.makedirs(MSH_DIR, exist_ok=True)

KEEP_SURFACE_GROUPS = {"inlet", "outlet", "walls"}

def convert(med_path, msh_path):
    m = meshio.read(med_path)

    # Build old_key -> new_tag mapping from cell_tags dict
    # cell_tags: {int_key: [name1, name2, ...]}
    old_to_new = {}  # old int key -> new gmsh physical tag
    surface_names = {}  # new_tag -> name (for field_data)

    new_tag = 1  # tag 1 = volume (unnamed, gmshToFoam won't make a patch)
    for key, names in m.cell_tags.items():
        k = int(key)
        # Volume group
        if any(n in ("Group_Of_All_Volumes", "fluid") for n in names):
            old_to_new[k] = 1
        # Surface patches we want
        for name in names:
            if name in KEEP_SURFACE_GROUPS:
                if name not in surface_names.values():
                    new_tag += 1
                    old_to_new[k] = new_tag
                    surface_names[new_tag] = name

    print("  tag remap:", old_to_new)
    print("  surfaces: ", surface_names)

    # Remap cell_data arrays: unknown tags -> 0, known -> new tag
    new_cell_data = []
    for arr in m.cell_data["cell_tags"]:
        remapped = np.zeros_like(arr)
        for old_key, new_t in old_to_new.items():
            remapped[arr == old_key] = new_t
        new_cell_data.append(remapped)

    m.cell_data["gmsh:physical"] = new_cell_data
    del m.cell_data["cell_tags"]

    # field_data tells gmsh writer the name+dim of each physical tag
    # Volume tag 1 is intentionally absent (unnamed = no patch in OpenFOAM)
    m.field_data = {
        name: np.array([tag, 2])   # dim=2 for surface patches
        for tag, name in surface_names.items()
    }

    meshio.write(msh_path, m, file_format="gmsh22", binary=False)
    print("[OK] {}".format(os.path.basename(msh_path)))


med_files = sorted(f for f in os.listdir(MED_DIR) if f.endswith(".med"))
if not med_files:
    print("No .med files found in: {}".format(MED_DIR))
    sys.exit(0)

print("Converting {} file(s)...\n".format(len(med_files)))
for fname in med_files:
    try:
        convert(os.path.join(MED_DIR, fname),
                os.path.join(MSH_DIR, fname.replace(".med", ".msh")))
    except Exception as e:
        import traceback
        print("[ERROR] {}: {}".format(fname, e))
        traceback.print_exc()

print("\nDone. MSH -> {}".format(MSH_DIR))