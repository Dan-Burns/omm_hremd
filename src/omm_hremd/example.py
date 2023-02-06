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
from omm_hremd.replex_tools import *

base_dir = '/Users/danielburns/Library/CloudStorage/Box-Box'
output_dir = '/Users/danielburns/Library/CloudStorage/Box-Box/notebooks/output_test'
structure = f'{base_dir}/notebooks/omm_replex_sandbox/output/1zip_10ns.pdb'
system = f'{base_dir}/notebooks/omm_replex_sandbox/output/1zip_system.xml'
n_replicas = 28
base_temperature = 280
lambdas = generate_lambdas(base_temperature,380,n_replicas)

# make an initial replex and equilibrate it with a small 
# timestep then use the sampler states to start a new one
# with a normal timestep
hremd_eq = make_hremd_simulation(structure,
                        system,
                        temperature=base_temperature*kelvin,
                        friction=100/picosecond,
                        output=f'{output_dir}/min_eq.nc',
                        timestep=0.01*femtosecond,
                        lambda_values = lambdas
                        )

hremd_eq.minimize()
hremd_eq.equilibrate(5)

equilibrated_states = hremd.sampler_states
hremd =  make_hremd_simulation(structure,
                        system,
                        temperature=base_temperature*kelvin,
                        friction=10/picosecond,
                        output=f'{output_dir}/run_output.nc',
                        timestep=2*femtosecond,
                        lambda_values = lambdas
                        )