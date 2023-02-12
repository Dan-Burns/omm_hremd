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
from cg_openmm.simulation.rep_exch import extract_trajectory

import shutil

run = 1
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


## Need to change rep_exch.py to allow for saving trajectories to a specified directory
## for now can just shutile.move() to specified directory
## I guess make_replica_dcd_files are the continuous trajectories that move up and down through temperature..
make_state_dcd_files(
    topology, timestep=timestep, time_interval=5,
    output_dir=output_dir, output_data="run_output.nc", checkpoint_data="run_output_checkpoint.nc",
    frame_begin=0, frame_stride=1, center=False):

# trajectory files will be labeled f"{output_dir}/replica_{replica_index+1}.dcd"

# Remove Waters and Center
xtc_directory = 'processed_xtc_files'
# place to store the xtc trajectories
make_directory(xtc_directory)

# move the generically named dcd files that cg_openmm.rep_exch.make_replica_dcd_files produced to a separate dir and add the run number to them
old_dcds = 'old_dcd_files'
make_directory(old_dcds)

for replica_index in range(n_replicas):
    dcd_file = f"{output_dir}/replica_{replica_index+1}.dcd"
    center_traj(structure, dcd_file, output=f'{xtc_directory}/traj_{run}_rep_{replica_index+1}.xtc', remove_waters=True,
                selection_string='protein')
    shutil.move(dcd_file, f'{old_dcds}/run_{run}_replica_{replica_index+1}.dcd')




    