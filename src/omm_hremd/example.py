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
from cg_openmm.rep_exch import extract_trajectory

base_dir = '/Users/danielburns/Library/CloudStorage/Box-Box'
output_dir = '/Users/danielburns/Library/CloudStorage/Box-Box/notebooks/output_test'
structure = f'{base_dir}/notebooks/omm_replex_sandbox/output/1zip_10ns.pdb'
system = f'{base_dir}/notebooks/omm_replex_sandbox/output/1zip_system.xml'
n_replicas = 28
base_temperature = 280
lambdas = generate_lambdas(base_temperature,380,n_replicas)
timestep = 2*femtosecond
output_interval = 5 # 5 exchange attempts X 1000 steps per attempt = 10 ns

topology, positions = PDBFile(structure).topology, PDBFile(structure).positions

## Save Gromacs files
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

equilibrated_states = hremd_eq.sampler_states
hremd =  make_hremd_simulation(structure,
                        system,
                        temperature=base_temperature*kelvin,
                        friction=10/picosecond,
                        output=f'{output_dir}/run_output.nc',
                        timestep=timestep,
                        lambda_values = lambdas
                        )

hremd.run(1000000)


make_replica_dcd_files(
    topology, timestep=timestep, time_interval=5,
    output_dir=output_dir, output_data="run_output.nc", checkpoint_data="run_output_checkpoint.nc",
    frame_begin=0, frame_stride=1, center=False):

# trajectory files will be labeled f"{output_dir}/replica_{replica_index+1}.dcd"

# Remove Waters and Center

    