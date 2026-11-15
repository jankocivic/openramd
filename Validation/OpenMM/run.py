from openmm import *
from openmm.app import *
from openmm.unit import *
from ramd import RAMD

# System setup
pdb = PDBFile("../system.pdb")
prmtop = AmberPrmtopFile("../system.parm7", periodicBoxVectors=pdb.topology.getPeriodicBoxVectors())

topology = prmtop.topology
positions = pdb.positions

# System Configuration
nonbondedMethod = PME
nonbondedCutoff = 0.9*nanometers
constraints = HBonds
rigidWater = True
constraintTolerance = 0.000001
ewaldErrorTolerance = 0.0005

# Integration Options
dt = 0.002*picoseconds
temperature = 300
friction = 1.0/picosecond
pressure = 1*bar
total_steps = 10000000

# Build system
system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff, constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
integrator = LangevinMiddleIntegrator(temperature*kelvin, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
barostat = system.addForce(MonteCarloBarostat(pressure, temperature))

# Build simulation
simulation = Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature*kelvin)

# Append reporters
dataReporter = StateDataReporter('log.txt', 10000, totalSteps=total_steps, step=True, speed=True, progress=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True)
dcdReporter = DCDReporter('trajectory.dcd', 10000)
simulation.reporters.append(dataReporter)
simulation.reporters.append(dcdReporter)

# Run
ligand = [48601, 48602]
receptor = [973]

ramd = RAMD(simulation, ligand=ligand, receptor=receptor,
            magnitude=2 * kilocalorie_per_mole / angstrom,
            log_file="ramd.log", r_min=0.1 * angstrom, r_max=30 * angstrom, seed=None)
exit_step = ramd.run(max_steps=total_steps)
print("exit time:", exit_step * 2 * femtosecond)
