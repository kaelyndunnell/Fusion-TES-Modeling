import os
import shutil


def clean_number_folders(directory_path):
    all_numbers = []

    # scan directory and find all strictly numeric folder names
    for entry in os.scandir(directory_path):
        if entry.is_dir():
            # .isdigit() checks if the string contains only digits
            if entry.name.isdigit():
                all_numbers.append(int(entry.name))

    if not all_numbers:
        print("No numbered folders found.")
        return

    # determine which folders to keep; always keep 0 and maximum num in list
    max_num = max(all_numbers)
    dirs_to_keep = {0, max_num}

    print(f"Keeping folder '0' and folder '{max_num}' (the maximum).")

    # delete the rest
    for entry in os.scandir(directory_path):
        if entry.is_dir() and entry.name.isdigit():
            folder_value = int(entry.name)

            if folder_value not in dirs_to_keep:
                folder_path = entry.path
                try:
                    # shutil.rmtree() deletes the directory and all its contents
                    shutil.rmtree(folder_path)
                    print(f"Deleted folder: {folder_path}")
                except Exception as e:
                    print(f"Error deleting {folder_path}: {e}")


if __name__ == "__main__":
    folders = [
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.32m_s_r0.13_l0.51",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.33m_s_r0.10_l0.75",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.33m_s_r0.13_l0.74",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.34m_s_r0.05_l0.72",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.35m_s_r0.10_l0.91",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.35m_s_r0.12_l1.09",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.35m_s_r0.20_l0.77",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.36m_s_r0.04_l0.68",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.36m_s_r0.09_l0.69",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.37m_s_r0.11_l0.99",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.39m_s_r0.05_l0.79",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.39m_s_r0.06_l0.79",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.40m_s_r0.02_l0.53",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.40m_s_r0.12_l1.06",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.41m_s_r0.16_l1.03",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.42m_s_r0.07_l0.89",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.42m_s_r0.14_l0.77",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.42m_s_r0.15_l0.71",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.43m_s_r0.15_l0.97",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.44m_s_r0.18_l0.93",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.44m_s_r0.19_l0.93",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.44m_s_r0.19_l1.01",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.45m_s_r0.05_l0.54",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.45m_s_r0.07_l0.53",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.45m_s_r0.16_l0.67",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.45m_s_r0.19_l0.96",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.46m_s_r0.08_l1.01",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.46m_s_r0.20_l0.98",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.47m_s_r0.03_l0.63",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.47m_s_r0.17_l0.94",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.47m_s_r0.17_l1.03",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.48m_s_r0.06_l0.68",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.48m_s_r0.09_l0.87",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.48m_s_r0.10_l0.56",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.48m_s_r0.11_l1.05",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.48m_s_r0.18_l0.73",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.49m_s_r0.05_l0.86",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.49m_s_r0.08_l1.07",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.49m_s_r0.14_l0.75",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.50m_s_r0.17_l1.00",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.50m_s_r0.19_l0.87",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.51m_s_r0.15_l0.94",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.51m_s_r0.16_l0.74",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.51m_s_r0.17_l0.60",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.53m_s_r0.08_l0.82",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.53m_s_r0.08_l0.97",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.53m_s_r0.19_l0.74",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.54m_s_r0.07_l0.57",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.54m_s_r0.09_l0.53",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.54m_s_r0.09_l0.94",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.54m_s_r0.16_l0.70",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.54m_s_r0.17_l0.64",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.55m_s_r0.08_l0.56",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.55m_s_r0.08_l0.64",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.55m_s_r0.17_l0.67",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.56m_s_r0.06_l0.56",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.56m_s_r0.18_l0.83",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.57m_s_r0.04_l1.08",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.57m_s_r0.07_l0.82",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.57m_s_r0.08_l1.01",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.57m_s_r0.10_l0.84",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.57m_s_r0.14_l0.61",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.58m_s_r0.03_l1.09",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.59m_s_r0.06_l0.76",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.60m_s_r0.03_l0.98",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.60m_s_r0.08_l0.94",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.61m_s_r0.06_l0.68",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.61m_s_r0.06_l0.79",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.62m_s_r0.07_l0.56",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.62m_s_r0.09_l0.86",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.63m_s_r0.08_l1.07",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.63m_s_r0.09_l0.88",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.63m_s_r0.15_l0.55",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.63m_s_r0.15_l0.83",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.64m_s_r0.11_l0.54",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.64m_s_r0.20_l0.57",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.65m_s_r0.06_l0.62",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.65m_s_r0.09_l0.54",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.66m_s_r0.08_l0.89",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.66m_s_r0.19_l1.02",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.67m_s_r0.06_l0.79",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.67m_s_r0.10_l0.65",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.67m_s_r0.15_l0.73",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.68m_s_r0.03_l1.09",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.68m_s_r0.04_l1.00",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.68m_s_r0.05_l0.61",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.69m_s_r0.11_l0.52",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.69m_s_r0.13_l0.59",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.69m_s_r0.17_l0.73",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.69m_s_r0.18_l1.10",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.69m_s_r0.20_l0.95",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.70m_s_r0.06_l0.75",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.70m_s_r0.12_l0.71",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.71m_s_r0.13_l0.85",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.72m_s_r0.09_l0.58",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.72m_s_r0.09_l0.64",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.73m_s_r0.08_l0.99",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.73m_s_r0.10_l0.85",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.73m_s_r0.20_l0.72",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.74m_s_r0.05_l0.96",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.75m_s_r0.04_l0.73",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.75m_s_r0.05_l0.90",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.75m_s_r0.11_l0.76",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.76m_s_r0.14_l0.82",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.76m_s_r0.19_l0.68",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.77m_s_r0.10_l0.80",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.77m_s_r0.11_l0.86",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.77m_s_r0.13_l0.83",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.78m_s_r0.06_l0.82",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.78m_s_r0.13_l0.99",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.78m_s_r0.19_l0.69",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.79m_s_r0.03_l0.58",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.80m_s_r0.03_l1.06",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.80m_s_r0.14_l0.58",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.81m_s_r0.05_l1.06",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.81m_s_r0.10_l0.65",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.81m_s_r0.14_l0.98",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.82m_s_r0.03_l0.73",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.82m_s_r0.18_l0.96",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.82m_s_r0.19_l0.59",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.83m_s_r0.13_l0.53",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.84m_s_r0.04_l0.98",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.84m_s_r0.11_l0.51",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.84m_s_r0.15_l1.09",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.84m_s_r0.18_l0.54",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.85m_s_r0.18_l1.10",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.85m_s_r0.19_l0.79",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.86m_s_r0.05_l0.65",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.86m_s_r0.09_l1.00",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.86m_s_r0.13_l0.95",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.86m_s_r0.17_l1.01",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.87m_s_r0.19_l0.97",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.88m_s_r0.03_l0.73",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.88m_s_r0.07_l0.76",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.88m_s_r0.09_l1.08",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.88m_s_r0.17_l0.94",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.89m_s_r0.02_l0.83",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.89m_s_r0.03_l1.03",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.89m_s_r0.04_l1.09",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.90m_s_r0.09_l0.55",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.90m_s_r0.11_l1.07",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.90m_s_r0.13_l0.88",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.90m_s_r0.14_l0.84",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.91m_s_r0.12_l0.87",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.91m_s_r0.17_l0.81",
        "/home/kaelyn/Fusion-TES-Modeling/openfoam/velocity_parametrization/lipb/pipe_0.92m_s_r0.04_l0.67",
    ]
    for TARGET_DIRECTORY in folders:
        clean_number_folders(TARGET_DIRECTORY)
