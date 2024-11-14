import numpy as np
import os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

# Define the van der Waals radii for relevant atoms (in angstroms)
vdw_radii = {
    'C': 1.70,  # Carbon
    'N': 1.55,  # Nitrogen
    'O': 1.52,  # Oxygen
    'H': 1.20,  # Hydrogen
    'F': 1.47,  # Fluorine
    'Cl': 1.75, # Chlorine
    'Br': 1.85, # Bromine
    'I': 1.98   # Iodine
}

# Fixed radius for Ba atoms (ghost atoms)
fixed_radius = 2.0

def read_xyz(file_path):
    """Read an XYZ file and return a list of atom data."""
    atoms = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines[2:]:  # Skip the first two lines (header)
            parts = line.split()
            atom_type = parts[0]
            x, y, z = map(float, parts[1:])
            atoms.append((atom_type, x, y, z))
    return atoms

def calculate_void_volume(ba_position, atoms, vdw_radii, grid_resolution=0.1):
    """Calculate the void volume around the Ba atom using a grid-based approach."""
    ba_radius = fixed_radius
    total_volume = 4/3 * np.pi * ba_radius**3  # Volume of sphere around Ba
    
    # Calculate the bounding box for the grid
    min_x = ba_position[0] - ba_radius
    max_x = ba_position[0] + ba_radius
    min_y = ba_position[1] - ba_radius
    max_y = ba_position[1] + ba_radius
    min_z = ba_position[2] - ba_radius
    max_z = ba_position[2] + ba_radius

    # Define the grid
    x_points = np.arange(min_x, max_x, grid_resolution)
    y_points = np.arange(min_y, max_y, grid_resolution)
    z_points = np.arange(min_z, max_z, grid_resolution)

    # Initialize a grid to count points
    occupied_points = 0
    total_points = 0

    for x in x_points:
        for y in y_points:
            for z in z_points:
                # For each grid point, check if it lies within the Ba sphere
                distance_to_ba = np.linalg.norm(np.array([ba_position[0] - x, ba_position[1] - y, ba_position[2] - z]))
                if distance_to_ba <= ba_radius:
                    total_points += 1
                    # Check if the point is inside any other atoms' radius
                    is_occupied = False
                    for atom_type, ax, ay, az in atoms:
                        atom_radius = vdw_radii.get(atom_type, 1.20)  # Default to H if not in dictionary
                        distance_to_atom = np.linalg.norm(np.array([ax - x, ay - y, az - z]))
                        if distance_to_atom <= atom_radius:
                            is_occupied = True
                            break  # No need to check further atoms if the point is occupied
                    if is_occupied:
                        occupied_points += 1

    # Calculate the void volume based on unoccupied grid points
    void_volume = total_volume * (1 - (occupied_points / total_points)) if total_points > 0 else 0
    return void_volume

def process_single_file(file_path):
    """Process a single XYZ file and calculate the void volumes for each Ba atom."""
    atoms = read_xyz(file_path)
    void_volumes = []
    
    # Identify Ba atoms and calculate the void volume around each Ba
    for atom_type, x, y, z in atoms:
        if atom_type == 'Ba':  # If the atom is Ba, treat it as a ghost atom
            ba_position = (x, y, z)
            void_volume = calculate_void_volume(ba_position, atoms, vdw_radii)
            void_volumes.append((file_path, void_volume))
    
    return void_volumes

def process_xyz_files_multithreaded(directory):
    """Process all XYZ files in the given directory using multiple threads."""
    void_volumes = defaultdict(list)
    
    with ThreadPoolExecutor() as executor:
        futures = []
        for filename in os.listdir(directory):
            if filename.endswith('.xyz'):
                file_path = os.path.join(directory, filename)
                futures.append(executor.submit(process_single_file, file_path))
        
        for future in futures:
            result = future.result()
            for filename, void_volume in result:
                # Extract the part of the file name after the '-' and group by that
                group_key = filename.split('-')[-1].split('.')[0]  # Use the part after '-' before '.xyz'
                void_volumes[group_key].append((filename, void_volume))

    return void_volumes

# Directory containing the .xyz files
directory = './pedaling_ghost/'

# Calculate void volumes and group by file name suffix using multithreading
void_volumes = process_xyz_files_multithreaded(directory)

# Sort and print the void volumes in ascending order
for group_key, volumes in void_volumes.items():
    # Sort by void volume
    sorted_volumes = sorted(volumes, key=lambda x: x[1])  # Sort by void volume (second element)
    print(f"Group: {group_key}")
    for filename, void_volume in sorted_volumes:
        print(f"  {filename}: {void_volume:.4f}")
