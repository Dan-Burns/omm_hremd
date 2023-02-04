import openmm
from openmm import *
from openmm.app import *
from openmm.unit import *
import parmed
import numpy   as np
import pandas  as pd
import MDAnalysis as mda
import openmmtools
from openmmtools import *
from openmmtools.multistate import MultiStateReporter, MultiStateSampler, ReplicaExchangeAnalyzer, ReplicaExchangeSampler
from sys import stdout



##### Helper Functions ########
#### Can move these to a utils.py ######
# generate the lambdas
def geometric_progression(tmin, tmax, n):
    '''
    make a geometric progression of temperatures in attempt to keep energy overlap consistent
    between all replicas and consequently equal exchange rates

    tmin: int or float
        Lowest temperature that replicas will be run at.

    tmax: int or float
        Maximum temperature that replicas will be run at.

    n: int
        Number of replicas in the replica exchange simulation. Should be a multiple of the number
        of available nodes.

    '''
    t = lambda tmin, tmax, n, i: tmin*np.exp(i*np.log(tmax/tmin)/(n-1))

    temps = []

    for i in range(n):
        temps.append(t(tmin,tmax,n,i))

    return temps

def generate_lambdas(tmin, tmax, n):
    '''
    Turn a temperature list into a set of lambdas for HREMD
    '''
    temps = geometric_progression(300, 320, n)
    return [min(temps)/temp for temp in temps]

class HREMDSimulation():
    '''
    Base class for constructing Hamiltonian Replica Exchange Simulation object.
    
    Parameters
    ----------
    structure : PDB file
        The PDB file written should be a pre-equilibrated system written from an openmm simulation.

    system : .xml file
        The system containing composed of a forcefield, nonbonded method, nonbonded cutoff,
        constraints, barostat, and temperature definitions
    
    
    '''

    def __init__(self, structure=None, system=None, 
                    box_vectors=None,
                    temperature=300 * kelvin,
                    integrator='LangevinSplittingDynamicsMove',
                    friction=10/picosecond,
                    platform='None',
                    output='output.nc',
                    output_interval=1,
                    exchange_attempts=1000000000,
                    steps_per_move=1000,
                    step_size=2*femtosecond,
                    lambda_values=None,
                    effective_temperatures=None,
                    lambda_selection='protein',
                    lambda_terms=['electrostatics,sterics,torsions'],
                    run_time=200*nanosecond,
                    ):



    def make_alchemical_system(self):
        '''
        Return a ReplicaExchangeSimulation object
        '''
        pdb = PDBFile(self.structure)
        positions, topology = pdb.positions, pdb.topology

        # open the system
        with open(self.system) as infile:
            system = XmlSerializer.deserialize(infile.read())

        ### Parmed uses amber masks for selection strings or you can turn the system into a df
        ### pmd_sys = parmed.openmm.load_topology(pdb.topology, system=system, xyz=pdb.positions)
        ### sys_df = pmd_sys.to_dataframe()     
        ### and select with something like list(sys_df[sys_df['chain']=='A'].index)
        ### I'd prefer using mdanalysis to select the alchemical atoms

        ### Get the atom indices that will have the reduced potential applied
        u = mda.Universe(pdb)
        selection = u.select_atoms(lambda_selection)
        # np.array of 0 indexed atom indices
        alchemical_atoms = selection.atoms.ix

        # if applying lambda to torsions then AlchemicalRegion(alchemical_torsions=True)
        alchemical_region = alchemy.AlchemicalRegion(alchemical_atoms=guest_atoms,alchemical_torsions=True)
        factory           = alchemy.AbsoluteAlchemicalFactory()
        alchemical_system = factory.create_alchemical_system(system, alchemical_region)
        alchemical_state  = alchemy.AlchemicalState.from_system(alchemical_system)

        # Initialize compound thermodynamic states at different temperatures and alchemical states.
        # for the Hamiltonian approach - same temperature but different lambdas.
        protocol = {'temperature': [300 for i in range(n_replicas)] * kelvin,
                    'lambda_electrostatics': lambdas,
                    'lambda_sterics': lambdas,
                    'lambda_torsions': lambdas}

        thermodynamic_states = states.create_thermodynamic_state_protocol(alchemical_system,
                                                                        protocol=protocol,
                                                                        composable_states=[alchemical_state])

        box_vectors = system.getDefaultPeriodicBoxVectors()

        # Generate the sampler states that will have their positions updated throughought the simulation
        sampler_states = [states.SamplerState(positions=positions, box_vectors=box_vectors) 
                  for i,_ in enumerate(thermodynamic_states)]
        
        
        # This is the simulation block that occurs between each swap attempt
        # if the simulation errors out, drop timestep to ~ 0.01 femtoseconds
        # and increase collision_rate to ~ 100/picosecond
        # and then run for ~100 moves to get the system to calm down
        move = mcmc.LangevinSplittingDynamicsMove(timestep=2*femtosecond,
                                                        n_steps=steps_per_move,
                                                collision_rate=friction
                                                        )

    

        simulation = ReplicaExchangeSampler(
                mcmc_moves=move,
                number_of_iterations=exchange_attempts,
                replica_mixing_scheme='swap-neighbors',
                )

        reporter = MultiStateReporter(output, checkpoint_interval=1)
        
        simulation.create(thermodynamic_states, sampler_states, reporter)

        

    def minimize_equilibrate():
        '''
        These HREMD sims seem to be really sensitive to initial conditions.
        An equilibrated input file will still blow up when the lambda is applied so this 
        is intended to try different equilibration timesteps and collision rates to get the systems happy before 
        moving it to run status.

        Equilibration will need to be confirmed with some metric (potential energy distribution) so that 
        test_exchange_rate() results can be considered reliable.
        '''

    def test_exchange_rate(self, target_exchange_rate=0.10)
        '''
        do multiple short runs and gradient descent towards a target exchange rate
        d_exchange_rate wrt number of systems
        repeatedly runs make_multi_sim() runs short mds
        target_exchange_rate of 0.10 by default because the goal is to sample the effective temperature space
        while minimizing computational cost as opposed to sampling conformations at one temperature which would benefit from higher
        exchange rates.
        '''

    def make_gromacs_files():
        '''
        use for post processing and recentering trajectories with gmx trjconv
        '''

class ReplicaExchangeOutputAnalyzer():
    '''
    Pull all the data/trajectories out of the Replica exchange .nc files.  Deal with multiple output files from restarts, organize them,
    rewrap/center the trajectories.  Read log files, plot exchange rates and track systems through lambda space.

    ''' 