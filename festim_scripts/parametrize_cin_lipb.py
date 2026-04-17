from main import build_festim_model

c_out_all = {}

c_in_atms_m3 = [
    1e15,
    1e16,
    1e17,
    1e18,
    1e19,
    1e20,
    1e21,
    1e22,
    1e23,
    1e24,
    1e25,
]  # atms/m3

for conc in c_in_atms_m3:

    print(f"Running simulation for inlet concentration of {conc} #/m3 in FliBe.")

    # LIPB
    my_model, c_out = build_festim_model(
        c_inlet=conc,
        residual_pressure=0,
        breeder="lipb",
        openfoam_data_file="openfoam/lipb_simple/tes.foam",  # TODO: CHANGE THIS
        openfoam_final_time=5600,
        festim_mesh_file="meshing/test_festim.msh",
        results_folder=f"lipb_festim_results/conc_{conc}",
    )

    # INITIALISE AND RUN
    my_model.initialise()
    my_model.run()

    c_out_all[conc] = c_out.data[0]

    print(f"Simulation for inlet concentration {conc} #/m3 complete.")
