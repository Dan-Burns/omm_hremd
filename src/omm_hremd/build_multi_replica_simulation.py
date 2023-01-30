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

class HREMDSimulation():
    '''
    Base class for constructing Hamiltonian Replica Exchange Simulation object.
    
    
    
    '''

    def __init__(self, pdb_file=None, system=None, 
                    box_vectors=None,
                    temperature=300 * kelvin,
                    integrator='LangevinSplittingDynamicsMove',
                    platform='None',
                    output='output.nc',
                    output_interval=1,
                    steps_per_move=1000,
                    step_size=2*femtosecond,
                    lambda_values=None,
                    effective_temperatures=None,
                    lambda_atoms='protein',
                    lambda_terms=['electrostatics,sterics,torsions'],
                    run_time=200*nanosecond,
                    target_exchange_rate=0.10)



    def make_multi_sim():
        '''
        Return a ReplicaExchangeSimulation object
        '''

    def minimize_equilibrate():
        '''
        These HREMD sims seem to be really sensitive to initial conditions.
        An equilibrated input file will still blow up when the lambda is applied so this 
        is intended to try different equilibration timesteps and collision rates to get the systems happy before 
        moving it to run status.

        Equilibration will need to be confirmed with some metric (potential energy distribution) so that 
        test_exchange_rate() results can be considered reliable.
        '''

    def test_exchange_rate(target_exchange_rate)
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