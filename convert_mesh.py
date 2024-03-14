import meshio


def convert_med_to_xdmf(filename, cell_type="tetra", facet_type="triangle"):
    """_summary_

    Args:
        filename (str): the filename of the .med file (must end with .med)
        cell_type (str, optional): the type of cell. Defaults to "tetra".
        facet_type (str, optional): the type of facet. Defaults to "triangle".

    Returns:
        dict: a correspondance dict from tag number to subdomain
    """
    
    mesh = meshio.read(filename)

    for mesh_block in mesh.cells:
        if mesh_block.type == cell_type:
            print(mesh.cell_data_dict["cell_tags"][cell_type])
            meshio.write_points_cells(
                "mesh_cells.xdmf",
                mesh.points,
                [mesh_block],
                cell_data={"f": [-1 * mesh.cell_data_dict["cell_tags"][cell_type]]},
            )
        elif mesh_block.type == facet_type:
            meshio.write_points_cells(
                "mesh_facets.xdmf",
                mesh.points,
                [mesh_block],
                cell_data={"f": [-1 * mesh.cell_data_dict["cell_tags"][facet_type]]},
            )

    correspondance_dict = mesh.cell_tags

    print(correspondance_dict)
    return correspondance_dict

if __name__ == "__main__":
    correspondance_dict = convert_med_to_xdmf("Mesh_1.med")
    print(correspondance_dict)