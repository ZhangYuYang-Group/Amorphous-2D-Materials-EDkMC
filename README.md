# EDkMC simulation of monolayer amorphous materials
A Energy-Driven kinetic monte-carlo algorithm for generating configurations of monolayer amorphous materials

## Requirements
To run this script, you need:
-python3 (numpy,pymatgen,scipy)
-LAMMPS (if you need LAMMPS calls PYTHON version, you need to compile the LAMMPS with PYTHON)

## How to use
Edit the "kMC_maBN.py" or "in.lammps" for your potential or specific parameters
then, simply run
PYTHON calls LAMMPS:
```
python kMC_maBN.py
```
LAMMPS calls PYTHON:
```
bash run.sh
```
or
```
lmp -i in.lammps
```

in a terminal to start the simulation

## Interatomic potentials

The empirical and machine-learning potentials for MAC and maBN are available in the [latest release](https://github.com/ZhangYuYang-Group/Amorphous-2D-Materials-EDkMC/releases/latest).

## How to cite
If you find this script userful, please cite:
```
Nano Lett. 2022, 22, 8018−8024
Phys. Rev. B 2024, 109, 174106
```
