#!/usr/bin/env python3
import sys
import torch
import torchani
import argparse
import numpy as np
from ase.io import read
from time import time

# Initialize Argument Parser
parser = argparse.ArgumentParser(description='Compute ligand internal energy.')
parser.add_argument('input_fn', type=str, help='Input file name')
parser.add_argument('output_fn', type=str, help='Output file name')
parser.add_argument('--use_gpu', action='store_true', help='Use GPU if available')

args = parser.parse_args()

# Determine device to use
device = torch.device('cuda' if torch.cuda.is_available() and args.use_gpu else 'cpu')

# Initialize ANI model
model = torchani.models.ANI2x(periodic_table_index=True).to(device)

# Conversion factor
hartree_to_kcal_per_mol = 627.5094740631

# Supported atomic numbers for ANI2x model
supported_atomic_numbers = {1, 6, 7, 8, 16, 9, 17}


def main():
    # Read the XYZ file into a list of Atoms objects
    molecules = read(args.input_fn, index=':')

    # Prepare output file
    with open(args.output_fn, 'w') as output:

        # Get species and coordinates
        species_list = [mol.get_atomic_numbers() for mol in molecules]
        coordinates_list = [mol.get_positions() for mol in molecules]

        # To fully utilize GPU, we could do two things:
        # The first way
        # Duplicate the list 10 times to get a 2000 batch size that could take advantage of GPU
        #    cpu: 1.2343s, gpu: 0.7428s
        # species_list = species_list * 10
        # coordinates_list = coordinates_list * 10

        # Check for unsupported elements
        if not all(elem in supported_atomic_numbers for mol_species in species_list for elem in mol_species):
            raise ValueError("Unsupported element detected")

        # Convert to tensors
        species = torch.tensor(np.array(species_list), device=device)
        coordinates = torch.tensor(np.array(coordinates_list), requires_grad=False, device=device, dtype=torch.float32)

        # Start timing
        if device.type == 'cuda':
            torch.cuda.synchronize()
        start_time = time()

        energy = hartree_to_kcal_per_mol * model((species, coordinates)).energies

        # To fully utilize GPU, we could do two things:
        # The second way
        # cuda initialization take a long time, so it is not fare to compare a single gpu run with cpu
        # Run 100 times you will see the difference, cpu 11.6s, gpu 2.3s
        # for i in range(100):
        #     # Compute energy
        #     energy = hartree_to_kcal_per_mol * model((species, coordinates)).energies

        # End timing
        if device.type == 'cuda':
            torch.cuda.synchronize()
        end_time = time()

        # Print result
        # Print result
        for en in energy:
            print(f'{en.item()}', file=output)

        print(f"Elapsed Time: {end_time - start_time:.4f} seconds")


if __name__ == '__main__':
    main()
