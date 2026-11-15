"""Random Acceleration Molecular Dynamics (RAMD) for OpenMM."""

import numpy as np
import openmm
from openmm import unit


class RAMD:
    """Drive an existing openmm.app.Simulation with a RAMD force."""

    def __init__(self, simulation, ligand, magnitude, receptor=None,
                 ramd_steps=50, r_min=0.025 * unit.angstrom,
                 r_max=50 * unit.angstrom, log_freq=50, log_file="ramd.log",
                 seed=None):
        self.sim = simulation
        self.ligand = list(ligand)
        self.receptor = list(receptor) if receptor else None
        self.ramd_steps = ramd_steps
        self.r_min = r_min.value_in_unit(unit.nanometer)
        self.r_max = r_max.value_in_unit(unit.nanometer)
        self.magnitude = magnitude.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        self.log_freq = log_freq
        self.log_file = log_file
        self.rng = np.random.default_rng(seed)
        self.direction = np.zeros(3)

        # Cache per-group atom indices and masses
        self._lig_idx = np.array(self.ligand)
        self._lig_mass = self._masses(self.ligand)
        self._rec_idx = np.array(self.receptor) if self.receptor else None
        self._rec_mass = self._masses(self.receptor) if self.receptor else None

        # A uniform force (fx, fy, fz) on the ligand centroid
        self.force = openmm.CustomCentroidBondForce(1, "-(fx*x1 + fy*y1 + fz*z1)")
        for name in ("fx", "fy", "fz"):
            self.force.addGlobalParameter(name, 0.0)
        self.force.addBond([self.force.addGroup(self.ligand)], [])
        self.sim.system.addForce(self.force)
        self.sim.context.reinitialize(preserveState=True)

    def _log(self, message):
        """Append one line to the log file."""
        with open(self.log_file, "a") as f:
            f.write(message + "\n")

    def _positions(self):
        return self.sim.context.getState(
            getPositions=True).getPositions(asNumpy=True)

    def _masses(self, indices):
        """Mass-weights (in daltons) for indices; computed once and cached."""
        return np.array(
            [self.sim.system.getParticleMass(i) / unit.dalton for i in indices]
        )

    @staticmethod
    def _com(positions, idx, masses):
        """Mass-weighted COM (nm) of idx from a positions Quantity."""
        coords = positions.value_in_unit(unit.nanometer)[idx]
        return np.average(coords, axis=0, weights=masses)

    def _direction_str(self):
        return np.array2string(self.direction, precision=3)

    def reorient(self):
        """Point the force in a fresh random direction at fixed magnitude."""
        self.direction = self.rng.normal(size=3)
        self.direction /= np.linalg.norm(self.direction)
        for name, value in zip(("fx", "fy", "fz"), self.magnitude * self.direction):
            self.sim.context.setParameter(name, value)

    def run(self, max_steps=100_000_000):
        """Run RAMD; return the step at which the ligand exited (or max_steps)."""
        open(self.log_file, "w").close()  # start with an empty log
        receptor_info = f"{len(self.receptor)} atoms" if self.receptor else "none"
        self._log(f"start: magnitude={self.magnitude:.3f} kJ/mol/nm "
                  f"ramd_steps={self.ramd_steps} r_min={self.r_min:.4f} nm "
                  f"r_max={self.r_max:.3f} nm ligand={len(self.ligand)} atoms "
                  f"receptor={receptor_info} openmm={openmm.version.full_version}")
        self.reorient()
        self._log(f"step = 0 initial direction = {self._direction_str()}")
        previous = self._com(self._positions(), self._lig_idx, self._lig_mass)
        steps = 0
        while steps < max_steps:
            self.sim.step(self.ramd_steps)
            steps += self.ramd_steps
            positions = self._positions()
            ligand_com = self._com(positions, self._lig_idx, self._lig_mass)
            displacement = np.linalg.norm(ligand_com - previous)
            distance = None
            if self.receptor is not None:
                receptor_com = self._com(positions, self._rec_idx, self._rec_mass)
                distance = np.linalg.norm(ligand_com - receptor_com)
            if steps % self.log_freq == 0:
                dist_str = f"{distance:.4f} nm" if distance is not None else "n/a"
                self._log(f"step = {steps} displacement = {displacement:.4f} nm "
                          f"distance = {dist_str}")
            if displacement < self.r_min:
                self.reorient()
                self._log(f"step = {steps} stalled ({displacement:.4f} < "
                          f"{self.r_min:.4f} nm): new direction = {self._direction_str()}")
            previous = ligand_com
            if distance is not None and distance > self.r_max:
                self._log(f"step = {steps} EXIT: ligand-receptor distance "
                          f"{distance:.4f} nm > {self.r_max:.4f} nm")
                return steps
        self._log(f"step = {steps} reached max_steps without an exit event")
        return steps
