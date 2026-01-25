conda install conda-forge::jupyter-book
conda install conda-forge::jupytext # necessary to open MyST files using BinderHub https://jupyterbook.org/en/stable/interactive/launchbuttons.html#launchbuttons-binder
conda install conda-forge::sphinx-tags
conda install conda-forge::matplotlib
conda install conda-forge::meshio[all]
conda install conda-forge::pip>=20.1
conda install conda-forge::ipykernel
conda install conda-forge::nbconvert
conda install conda-forge::fenics-dolfinx
conda install conda-forge::festim=2.0b1
conda install conda-forge::python-gmsh
conda install conda-forge::shapely
conda install conda-forge::scipy

pip install gmsh 
pip install pyparsing
pip install h-transport-materials==0.16
pip install pyvista[all,trame]  # trame cannot be installed on conda