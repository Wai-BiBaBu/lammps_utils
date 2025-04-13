"""
## LAMMPS Dump File Processor
This script processes all LAMMPS `*.dump` files (located in `INPUT_DIR`) from molecular dynamics simulations to extract specific timesteps and convert them to a single `.extxyz` format file (`OUTPUT_FILE`) with additional information including:
- Atomic positions
- Forces 
- Total energy (sum of individual atom energies)

## Input
All `*.dump` file located in `INPUT_DIR`.

## Output
One .extxyz file named `OUTPUT_FILE`.

## Why Not Use ASE?
While ASE (Atomic Simulation Environment) is a powerful toolkit, it's too slow when processing large datasets (e.g., 1,000,000 frames), even with 32 cores CPU. This custom script provides better performance for bulk processing.

## Frame Extraction Behavior
The script extracts frames at regular intervals:
- Default interval: every 100 timesteps (`TIMESTEP_INTERVAL`)
- Default range: from 0 up to 100,000 (`MAX_TIMESTEP`)

Example with default settings extracts frames:  
`0, 100, 200, ..., 100000`

## Input File Format Specification
The script expects LAMMPS dump files with the following ATOMS format:
ITEM: ATOMS id type x y z vx vy vz fx fy fz atom_energy

If not, line 110-114 should be changed.

## Atom keymap
`ATOM_TYPE_MAP` should be modified.

"""

import os
import glob
from multiprocessing import Pool

# Global configuration variables
INPUT_DIR = './lammps'              # Directory containing LAMMPS dump files
OUTPUT_FILE = 'all_per100.extxyz'   # Output file in extended XYZ format
MAX_TIMESTEP = 100000               # Maximum timestep to consider
TIMESTEP_INTERVAL = 100            # Interval between saved timesteps

# Atom type mapping dictionary
ATOM_TYPE_MAP = {
    "1": "C",    # Map type 1 to Carbon
    "2": "Si",   # Map type 2 to Silicon
    "3": "O"     # Map type 3 to Oxygen
    # Add more mappings as needed
}

def find_dump_files(root_dir):
    """Recursively find all .dump files in subdirectories of root_dir"""
    dump_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith('.dump'):
                dump_files.append(os.path.join(dirpath, filename))
    return dump_files

def process_dump_file(dump_file):
    """Process a single dump file, return content of all selected frames"""
    desired_timesteps = list(range(0, MAX_TIMESTEP + 1, TIMESTEP_INTERVAL))
    results = []
    
    with open(dump_file, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            
            if "ITEM: TIMESTEP" not in line:
                continue
                
            timestep = int(f.readline().strip())
            
            if timestep not in desired_timesteps:
                # Skip unwanted frames
                for _ in range(6): f.readline()  # Skip header sections
                line = f.readline()  # ITEM: ATOMS...
                if not line: break
                num_atoms = int(f.readline().split()[0])
                for _ in range(num_atoms):
                    f.readline()
                continue
                
            # Read number of atoms
            f.readline()  # ITEM: NUMBER OF ATOMS
            num_atoms = int(f.readline().strip())
            
            # Read box boundaries
            f.readline()  # ITEM: BOX BOUNDS pp pp pp
            xlo, xhi = map(float, f.readline().split())
            ylo, yhi = map(float, f.readline().split())
            zlo, zhi = map(float, f.readline().split())
            
            # Calculate box dimensions
            lx = xhi - xlo
            ly = yhi - ylo
            lz = zhi - zlo
            
            # Read atom data
            f.readline()  # ITEM: ATOMS ...
            atoms_data = []
            total_energy = 0.0
            for _ in range(num_atoms):
                parts = f.readline().split()
                atom_type = parts[1]  # First column is atom type
                mapped_type = ATOM_TYPE_MAP.get(atom_type, "X")  # Default to "X" if type not mapped
                x, y, z = map(float, parts[2:5])
                fx, fy, fz = map(float, parts[8:11])
                energy = float(parts[11])
                total_energy += energy
                atoms_data.append((mapped_type, x, y, z, fx, fy, fz))
            
            # Prepare output content
            base_name = os.path.splitext(os.path.basename(dump_file))[0]
            rel_path = os.path.relpath(os.path.dirname(dump_file), INPUT_DIR)
            category = f"{rel_path}/{base_name}-{timestep}" if rel_path != '.' else f"{base_name}-{timestep}"
            
            frame_content = [f"{num_atoms}"]
            frame_content.append(f'Lattice="{lx} 0.0 0.0 0.0 {ly} 0.0 0.0 0.0 {lz}" '
                              f'Properties=species:S:1:pos:R:3:forces:R:3 '
                              f'name=Amorphous_Bulk category={category} '
                              f'energy={total_energy} pbc="T T T"')
            
            for atom_type, x, y, z, fx, fy, fz in atoms_data:
                frame_content.append(
                    f"{atom_type} {x:.8f} {y:.8f} {z:.8f} {fx:.8f} {fy:.8f} {fz:.8f}"
                )
            
            results.append("\n".join(frame_content))
    
    return results

def main():
    # Recursively find all dump files
    dump_files = find_dump_files(INPUT_DIR)
    
    if not dump_files:
        print(f"No .dump files found in {INPUT_DIR} and its subdirectories")
        return
    
    print(f"Found {len(dump_files)} .dump files, starting processing...")
    print(f"Using atom type mapping: {ATOM_TYPE_MAP}")
    
    # Use multiprocessing for parallel processing
    with Pool() as pool:
        all_results = pool.map(process_dump_file, dump_files)
    
    # Write to unified output file
    with open(OUTPUT_FILE, 'w') as outfile:
        total_frames = 0
        for file_results in all_results:
            for frame in file_results:
                outfile.write(frame + "\n")
                total_frames += 1
    
    print(f"Processing complete! Processed {len(dump_files)} files, {total_frames} frames saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
