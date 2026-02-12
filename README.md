# Tritium Extraction Systems (TES) for Fusion Fuel Cycles

This repository develops a shell and tube tritium extractor design for extracting tritium from fusion liquid breeders. The shell and tube model is based on the design from [Conceptual design of a PAV-based tritium extractor for the WCLL breeding blanket of the EU DEMO: Effects of surface-limited vs. diffusion-limited modeling](https://doi.org/10.1016/j.fusengdes.2021.112363) by R. Bonifetto et al (2021). Using OpenFOAM for CFD simulations as coupled to FESTIM for tritium transport, the shell and tube extractor can be modeled and tested as a potential tritium extraction device. 

## How to Run: 

Clone the repository: 

 ```
git clone https://github.com/kaelyndunnell/Fusion-TES-Modeling
cd Fusion-TES-Modeling
```

This package requires two environments to run, one for CADQuery and one for OpenFOAM/FESTIM. To set up a new environment with the right dependencies for CADQuery, use: 

Then, set up a new environment with the right dependencies (e.g. dolfinx, FESTIM) using:

```
mamba create -f cadquery_env.yml 
```

Then, activate the environment to use CADQuery for geometry creation: 

```
mamba activate cadquery-env
```

Deactive the environment using 
```
conda deactivate
```

and set up a second environment with the right meshing, OpenFOAM (Foundation v13), and FESTIM dependencies with: 
```
conda env create -f environment.yml
```

Activate this environment using 
```
conda activate tes-pav-env
```
