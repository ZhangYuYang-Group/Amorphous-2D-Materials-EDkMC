# Energy-driven kinetic Monte Carlo simulation of monolayer amorphous materials

An energy-driven kinetic Monte Carlo (EDkMC) framework for generating and sampling atomic configurations of two-dimensional amorphous materials.

This repository provides the simulation codes associated with our studies of **monolayer amorphous carbon (MAC)** and **monolayer amorphous boron nitride (maBN)**. The method combines topology-changing Monte Carlo moves with structural relaxation and energy evaluation, enabling efficient exploration of the configurational space of large two-dimensional amorphous systems.

The original kMC framework and the resulting structural pictures of MAC and maBN were reported in our *Nano Letters* work. In our subsequent *Physical Review B* study, machine-learning potentials (MLPs) with density-functional-theory-level accuracy were introduced to generate MAC and maBN structures with different degrees of disorder and to investigate their structural and mechanical properties.

## Method overview

For each EDkMC step, the simulation generally follows the workflow:

1. Generate a trial configuration through a topology-changing Monte Carlo move.
2. Relax the trial structure and evaluate its energy using an interatomic potential.
3. Accept or reject the trial configuration according to the energy-based acceptance criterion.
4. Repeat the procedure to explore the configurational space and generate amorphous structures with different degrees of disorder.

Two implementations are provided in this repository:

* **Python-driven EDkMC:** Python controls the Monte Carlo procedure and calls LAMMPS for structural relaxation and energy evaluation.
* **LAMMPS-driven EDkMC:** LAMMPS performs the simulation and calls Python for the Monte Carlo operations.

EDkMC is used here primarily as a **configurational sampling and structural-search method**. The Monte Carlo step should not be interpreted as physical time, and the algorithmic acceptance parameter should not be directly identified with an experimental temperature or synthesis condition.

## Requirements

### Python

Python 3 with:

* NumPy
* SciPy
* pymatgen

### LAMMPS

LAMMPS is required for structural relaxation and energy evaluation.

For the **Python-driven** implementation, a working LAMMPS Python interface is required.

For the **LAMMPS-driven** implementation, LAMMPS should be compiled with the required Python support.

### DeePMD-kit

DeePMD-kit and its LAMMPS interface are required when using the machine-learning potentials `MAC.pb` and `maBN.pb`.

## How to use

Before running a simulation, edit `kMC_maBN.py` or `in.lammps` according to the material, interatomic potential, and EDkMC parameters of interest.

### Python-driven EDkMC

Run:

```bash
python3 kMC_maBN.py
```

### LAMMPS-driven EDkMC

Run the provided shell script:

```bash
bash run.sh
```

or directly run the LAMMPS input:

```bash
lmp -in in.lammps
```

The provided maBN example can be adapted to other monolayer amorphous systems by modifying the structure input, interatomic potential, trial moves, and corresponding simulation parameters.

## Interatomic potentials

The empirical and machine-learning potentials used in our MAC and maBN studies are available in the [latest release](https://github.com/ZhangYuYang-Group/Amorphous-2D-Materials-EDkMC/releases/latest).

| System | Potential                           | File        | Related work               |
| ------ | ----------------------------------- | ----------- | -------------------------- |
| MAC    | Empirical potential (AIREBO)        | `CH.airebo` | *Nano Letters* (2022)      |
| maBN   | Empirical potential                 | `BN.extep`  | *Nano Letters* (2022)      |
| MAC    | Machine-learning potential (DeePMD) | `MAC.pb`    | *Physical Review B* (2024) |
| maBN   | Machine-learning potential (DeePMD) | `maBN.pb`   | *Physical Review B* (2024) |

The empirical potentials were used in our original investigation of the contrasting atomic structures of elemental MAC and binary maBN. The DeePMD MLPs were subsequently developed to provide DFT-level accuracy at substantially reduced computational cost, allowing large-scale simulations of structures with different degrees of disorder and their mechanical responses.

For additional information on the potentials and their intended scope, see the description of the corresponding GitHub Release.

## References and citation

If you use this repository, please cite the corresponding work:

**[1] Structural framework and empirical-potential kMC simulations**

Y.-T. Zhang, Y.-P. Wang, X. Zhang, Y.-Y. Zhang, S. Du, and S. T. Pantelides,
“Structure of Amorphous Two-Dimensional Materials: Elemental Monolayer Amorphous Carbon versus Binary Monolayer Amorphous Boron Nitride,”
*Nano Letters* **22**, 8018–8024 (2022).
https://doi.org/10.1021/acs.nanolett.2c02542

**[2] Machine-learning-potential simulations and structural/mechanical properties**

X. Zhang, Y.-T. Zhang, Y.-P. Wang, S. Li, S. Du, Y.-Y. Zhang, and S. T. Pantelides,
“Structural and mechanical properties of monolayer amorphous carbon and boron nitride,”
*Physical Review B* **109**, 174106 (2024).
https://doi.org/10.1103/PhysRevB.109.174106

If the DeePMD potentials are used, please also cite DeePMD-kit as appropriate.
