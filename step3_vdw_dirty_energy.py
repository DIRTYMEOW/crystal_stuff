import os
import numpy as np

# Bondi's van der Waals radii for common atoms (in angstroms)
bondi_radii = {
    'C': 1.70,
    'H': 1.20,
    'N': 1.55,
    'O': 1.52,
    'F': 1.47,
    'Cl': 1.75,
    'Br': 1.85,
    'I': 1.98,
    'S': 1.80,
    'Ba': 2.68  # Add barium radius for reference
}

def read_xyz_file(filepath):
    """Read an XYZ file and return the atom types and their coordinates."""
    with open(filepath, 'r') as f:
        lines = f.readlines()[2:]  # Skip first two lines (header)
    coords = []
    for line in lines:
        parts = line.split()
        coords.append([float(parts[1]), float(parts[2]), float(parts[3]), parts[0]])
    return coords

def calculate_lj_potentials(coords, barium_coord):
    """Calculate the LJ potential energies contributed to the barium atom from all other atoms."""
    lj_potentials = []

    for coord in coords:
        element = coord[3]
        if element in bondi_radii:
            # Calculate radius based on Bondi's radii
            radius = bondi_radii[element] + bondi_radii['Ba']  # Using Ba radius as a reference
        else:
            print(f"What's that atom? Check if listed in Bondi, else add by yourself")
            continue

        # Calculate distance between the barium atom and the current atom
        dx = coord[0] - barium_coord[0]
        dy = coord[1] - barium_coord[1]
        dz = coord[2] - barium_coord[2]
        distance = (dx**2 + dy**2 + dz**2) ** 0.5

        # Calculate the LJ potential only if the distance is greater than zero
        if distance > 0:  # Ensure we don't divide by zero
            # Calculate the LJ potential, only using the (radius/distance)**12 term
            if distance < radius:
                lj_potential = (radius / distance) ** 12 - (radius / distance) ** 6
            else:
                lj_potential = 0.0

            lj_potentials.append(lj_potential)

    # Sum the LJ potentials
    lj_sum = sum(lj_potentials)
    return lj_sum

def main():
    total_energies = []  # List to hold total energies for each file
    # Process each .xyz file in the specified directory
    for filename in os.listdir('./pedaling_ghost/'):
        if filename.endswith('.xyz'):
            filepath = os.path.join('./pedaling_ghost/', filename)
            coords = read_xyz_file(filepath)

            # Identify the barium atom coordinates
            barium_coord = None
            for coord in coords:
                if coord[3] == 'Ba':
                    barium_coord = coord[:3]  # Extract coordinates (x, y, z)
                    break

            if barium_coord is None:
                print(f"No barium atom found in {filename}. Skipping.")
                continue

            # Calculate LJ potentials for the barium atom
            lj_sum = calculate_lj_potentials(coords, barium_coord)

            # Store the filename and its corresponding total energy
            total_energies.append((filename, lj_sum))

    # Sort total energies by the LJ potential values (ascending)
    total_energies.sort(key=lambda x: x[1])

    # Print the sorted total energies
    print("Total LJ potentials for Barium in ascending order:")
    for filename, lj_sum in total_energies:
        print(f"{filename}: {lj_sum:.6f}")

if __name__ == '__main__':
    main()
