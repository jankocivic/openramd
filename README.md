# openramd
This is a simple [OpenMM](https://openmm.org) plugin to run Random Acceleration Molecular Dynamics (RAMD).
It does not require installing any external dependencies.
The whole plugin is just one python file containing around 100 lines of code.
No noticeable overhead is introduced compared to running a simulation without RAMD.

RAMD is a technique that accelerates ligand unbinding from a binding site without biasing towards any specific direction.
If repeated many times it can be used to analyse ligand dissociation mechanisms and to compute the relative residence times of different ligands.
For more details on the technique and algorithm we suggest reading the original RAMD paper by [Lüdemann et al.](https://doi.org/10.1006/jmbi.2000.4154) and τRAMD paper by [Kokh et al.](https://doi.org/10.1021/acs.jctc.8b00230) which were used as a reference for this plugin.

Another similar python package already exists [openmm_ramd](https://github.com/seekrcentral/openmm_ramd) and it served as inspiration for this plugin, but here the whole package is condensed into one relatively short file without loss of any functionality.
In our opinion this makes it easier to maintain and use.

https://github.com/user-attachments/assets/3936fdce-4ca1-48ec-a141-3d3d5eda3300

## Installation

openramd needs [OpenMM](https://openmm.org), which it does not install for you.
Install OpenMM first (conda or pip, any version ≥ 7.6), then install openramd
into the same environment:

```bash
pip install git+https://github.com/jankocivic/openramd
```

or with uv:

```bash
uv pip install git+https://github.com/jankocivic/openramd   # into the active environment
uv add git+https://github.com/jankocivic/openramd           # into a uv project
```

Alternatively, skip installation and copy `ramd.py` next to your script.
Either way, import it with:

```python
from ramd import RAMD
```

## Usage

Set up and equilibrate your system as usual, then hand the `Simulation` to `RAMD`:

```python
ramd = RAMD(simulation, ligand=ligand, receptor=receptor,
            magnitude=14 * unit.kilocalorie_per_mole / unit.angstrom,
            log_file="ramd.log", seed=1)
exit_step = ramd.run()
print("exit time:", exit_step * 2 * unit.femtosecond)
```

`run()` returns the step at which the ligand exited. If it returns `max_steps`,
no exit occurred.

## Parameters

| Parameter    | Default        | Description                                                        |
|--------------|----------------|--------------------------------------------------------------------|
| `simulation` | —              | An `openmm.app.Simulation`. The RAMD force is added to its system. |
| `ligand`     | —              | Atom indices of the ligand (0-indexed).                            |
| `magnitude`  | —              | Force magnitude (OpenMM quantity, e.g. kcal/mol/Å).                |
| `receptor`   | `None`         | Atom indices of the receptor, used only to detect the exit.        |
| `ramd_steps` | `50`           | MD steps between direction checks.                                 |
| `r_min`      | `0.025 Å`      | Minimum ligand displacement per `ramd_steps` to keep the direction.|
| `r_max`      | `50 Å`         | Ligand–receptor COM distance that counts as an exit.               |
| `log_freq`   | `50`           | Logging interval in steps (use a multiple of `ramd_steps`).        |
| `log_file`   | `"ramd.log"`   | File the RAMD log is written to.                                   |
| `seed`       | `None`         | Seed for the random force directions.                              |

`magnitude`, `r_min` and `r_max` must be given as OpenMM quantities with units
(e.g. `0.1 * unit.angstrom`, `14 * unit.kilocalorie_per_mole / unit.angstrom`);
any compatible unit works and is converted internally. Plain numbers are not accepted.

`run(max_steps=100_000_000)` runs RAMD until exit or `max_steps`.

## Notes

- **Parameters.** The default parameters probably need to be adjusted for every
  system. Especially important are `magnitude`, `r_min` and `ramd_steps`.
- **Periodic boxes.** Distances are computed from OpenMM's unwrapped coordinates.
  Make sure the ligand and all receptor chains are whole and in the same periodic
  image in the starting structure, and choose `r_max` below half the shortest box length.
- **No receptor, no exit.** Without `receptor`, the run always continues to `max_steps`.
- **Potential energy.** The RAMD force adds a term to the reported potential
  energy that has no physical meaning.

## Validation

We compared openramd with the RAMD implementation in NAMD on oxygen unbinding
from the fluorescent protein KillerRed. With each code, we ran 100 RAMD
trajectories from the same starting structure with identical RAMD parameters:

| Parameter    | Value              |
|--------------|--------------------|
| `magnitude`  | 2 kcal mol⁻¹ Å⁻¹   |
| `r_min`      | 0.1 Å              |
| `ramd_steps` | 50                 |
| `r_max`      | 30 Å               |

The NAMD and OpenMM input files are in the `Validation` folder.

**Exit times.** Both codes give the same average exit-time with very similar distributions:

| N = 100              | NAMD       | openramd   |
|----------------------|------------|------------|
| Mean exit time (ns)  | 0.65  | 0.67  |
| Median exit time (ns)  | 0.54  | 0.44  |
| Standard deviation (ns)  | 0.49  | 0.64  |
| Minimum (ns)  | 0.06  | 0.06  |
| Maximum (ns)  | 2.82  | 3.34  |

A two-sample Kolmogorov–Smirnov test shows no significant difference between the two distributions (**p** = 0.37). The small differences could come from sampling noise or differences in the OpenMM and NAMD simulation engines not related with the RAMD implementation.

**Performance.** On a single NVIDIA L40 GPU and 16 core INTEL(R) XEON(R) GOLD 5515+, openramd runs at about 250 ns/day,
compared with about 100 ns/day for NAMD 3.0.2 on the same hardware: a
2.5× speedup (System: 63589 atoms, 2 fs, same trajectory writing frequency). We noticed that the GPU utilization
is 100% for OpenMM and 50% for NAMD. Possibly with an even more powerful CPU the performance of NAMD would approach OpenMM.
