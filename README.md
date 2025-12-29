# Tritium Extraction Systems (TES) for Fusion Fuel Cycles

This repository develops a shell and tube tritium extractor design for extracting tritium from fusion liquid breeders. The shell and tube model is based on the design from [Conceptual design of a PAV-based tritium extractor for the WCLL breeding blanket of the EU DEMO: Effects of surface-limited vs. diffusion-limited modeling](https://doi.org/10.1016/j.fusengdes.2021.112363) by R. Bonifetto et al (2021). Using OpenFOAM for CFD simulations as coupled to FESTIM for tritium transport, the shell and tube extractor can be modeled and tested as a potential tritium extraction device. 

## How to Run: 

Clone the repository: 

 ```
git clone https://github.com/kaelyndunnell/Fusion-TES-Modeling
cd Fusion-TES-Modeling
```

Run this command to create a new environment with the right dependencies:

```
conda env create -f environment.yml 
```

Then, activate the environment: 

```
conda activate tes-model-env
```