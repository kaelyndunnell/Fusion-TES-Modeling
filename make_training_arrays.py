import os
import re
import csv

FOLDER_RE = re.compile(
    r"in_(?P<in>[0-9eE+\-.]+)"
    r"_vel_(?P<vel>[0-9eE+\-.]+)"
    r"_(?:bend_|r)(?P<r>[0-9eE+\-.]+)"
    r"_(?:len_|l)(?P<l>[0-9eE+\-.]+)"
)

def last_value(filepath):
    last = None
    with open(filepath, newline="") as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            if len(row) >= 2 and row[1].strip():
                last = row[1].strip()
    if last is None:
        raise ValueError(f"No data rows found in {filepath}")
    return float(last)

def main():
    rows_x = []
    rows_y = []

    subfolders = sorted(os.listdir("lipb_festim_results/parametric_mesh"))
    for name in subfolders:
        folder_path = os.path.join("lipb_festim_results/parametric_mesh", name)
        m = FOLDER_RE.search(name)
        if not m:
            continue
        c_out_path = os.path.join(folder_path, "c_out.csv")
        flux_path  = os.path.join(folder_path, "permeation_flux.csv")
        c_out  = last_value(c_out_path)
        flux   = last_value(flux_path)
        rows_x.append([
            float(m.group("in")),
            float(m.group("vel")),
            float(m.group("r")),
            float(m.group("l")),
        ])
        rows_y.append([c_out, flux])

    if not rows_x:
        print("No valid simulation folders found. Nothing written.")
        return

    with open('lipb_emulators/parametric_mesh/inputs.csv', "w", newline="") as f:
        writer = csv.writer(f)
        # writer.writerow(["inlet_concentration", "velocity", "bend_radius", "length"])
        writer.writerows(rows_x)

    with open('lipb_emulators/parametric_mesh/outputs.csv', "w", newline="") as f:
        writer = csv.writer(f)
        # writer.writerow(["c_out", "permeation_flux"])
        writer.writerows(rows_y)

if __name__ == "__main__":
    main()