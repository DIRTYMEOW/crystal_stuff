import os
import numpy as np
from scipy.spatial import distance_matrix

# Directories
xyz_dir = './big_xyz'
output_dir = './pedaling_ghost'  # Directory to save the modified output files
carbon_nitrogen_threshold = 1.8       # Bond length threshold for C-N neighbors
halogen_bond_length_threshold = 2.2   # Extended bond length for halogens
sulfur_bond_length_threshold = 1.8    # Extended bond length for sulfur

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to calculate distance between two points
def calculate_distance(point1, point2):
    return np.linalg.norm(point1 - point2)

# Function to process each .xyz file and save the modified molecule
def process_xyz_file(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()
        
        # First line is the number of atoms, second line is a comment
        atom_count = int(lines[0].strip())
        atoms = []
        
        # Read atom coordinates
        for line in lines[2:2 + atom_count]:
            parts = line.split()
            atom_type = parts[0]
            x, y, z = map(float, parts[1:4])
            atoms.append((atom_type, np.array([x, y, z])))
        
        # Calculate the coordinates array and create a distance matrix
        coords_array = np.array([coord for _, coord in atoms])
        dist_matrix = distance_matrix(coords_array, coords_array)

        # Identify central C and N pairs within the specified distance threshold
        central_pairs = []
        for i, (atom1, coord1) in enumerate(atoms):
            if atom1 != 'C':
                continue
            for j, (atom2, coord2) in enumerate(atoms):
                if atom2 != 'N' or i == j:
                    continue
                distance = dist_matrix[i][j]
                if distance < carbon_nitrogen_threshold:
                    central_pairs.append((i, j))  # Store indices of central C and N

        if not central_pairs:
            print(f"No central C-N pair found in {filepath}")
            return

        # Choose the central C-N pair closest to the molecule's center
        center_coord = np.mean(coords_array, axis=0)
        central_pair = min(
            central_pairs,
            key=lambda pair: calculate_distance((atoms[pair[0]][1] + atoms[pair[1]][1]) / 2, center_coord)
        )
        central_c_index, central_n_index = central_pair

        # Use depth-first search to find the molecule connected to the central C and N
        visited = set()

        def find_connected_atoms(start_idx):
            """Recursively find all connected atoms within bonding thresholds."""
            molecule = []
            stack = [start_idx]
            while stack:
                idx = stack.pop()
                if idx not in visited:
                    visited.add(idx)
                    molecule.append(idx)
                    for neighbor_idx in range(len(atoms)):
                        if neighbor_idx == idx:
                            continue
                        distance = dist_matrix[idx][neighbor_idx]
                        # Determine threshold based on neighbor atom type
                        neighbor_atom_type = atoms[neighbor_idx][0]
                        if neighbor_atom_type in ['F', 'Cl', 'Br', 'I']:
                            threshold = halogen_bond_length_threshold
                        elif neighbor_atom_type == 'S':
                            threshold = sulfur_bond_length_threshold
                        else:
                            threshold = carbon_nitrogen_threshold
                        
                        if distance < threshold and neighbor_idx not in visited:
                            stack.append(neighbor_idx)
            return molecule

        # Get the full molecule connected to both central C and N atoms
        central_molecule = set(find_connected_atoms(central_c_index)).union(
            find_connected_atoms(central_n_index)
        )

        # Calculate the midpoint between the central C and N for the barium placement
        barium_position = (atoms[central_c_index][1] + atoms[central_n_index][1]) / 2

        # Create a new list excluding atoms in the central molecule
        modified_atoms = [(atom, coord) for i, (atom, coord) in enumerate(atoms) if i not in central_molecule]

        # Append the new barium atom at the midpoint
        modified_atoms.append(('Ba', barium_position))

        # Save the modified molecule to an .xyz file
        output_filepath = os.path.join(output_dir, f"{os.path.basename(filepath)}")
        with open(output_filepath, 'w') as output_file:
            output_file.write(f"{len(modified_atoms)}\n")
            output_file.write(f"Modified molecule with central C-N replaced by Ba from {os.path.basename(filepath)}\n")
            for atom_type, coord in modified_atoms:
                output_file.write(f"{atom_type} {coord[0]} {coord[1]} {coord[2]}\n")
        
        # Print confirmation and details
        print(f"Processed file: {os.path.basename(filepath)}")
        print(f"Modified molecule saved to {output_filepath}")

# Iterate through all .xyz files in the directory
for filename in os.listdir(xyz_dir):
    if filename.endswith('.xyz'):
        process_xyz_file(os.path.join(xyz_dir, filename))
