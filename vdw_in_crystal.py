# add or change whatever radii u want here, Bondi’s reference used here
bondi_radii = {
    'C': 1.70,
    'H': 1.20,
    'N': 1.55,
    'O': 1.52,
    'F': 1.47,
    'Cl': 1.75,
    'Br': 1.85,
    'I': 1.98,
    'S': 1.80
}

# pseudo atom coordinate is here~~~~~~~~~~remember to remove central molecule manually, and record the middle-point to C=N and input here
pseudo_coord = [0.0, 0.0, 0.0]
pseudo_radius = 2.5

# .xyz file name is here~~~~~~~~~~~~~~~ for other file, u need to re-write coding here
with open('*.xyz', 'r') as f:
    lines = f.readlines()[2:] # skip first two lines in xyz file to avoid iqmol output error
    coords = []
    for line in lines:
        parts = line.split()
        coords.append([float(parts[1]), float(parts[2]), float(parts[3]), parts[0]])

# loop for LJ potentials
lj_potentials = []
for coord in coords:
# calculate LJ sigma/epsilon
    element = coord[3]
    if element in bondi_radii:
        radius = bondi_radii[element] + pseudo_radius
        epsilon = (bondi_radii[element] * pseudo_radius) ** 0.5
    else:
        print(f"What's that atom? Check if listed in Bondi, else add by urself")
        continue

# calculate distance between pseudo-atom and current atom
    dx = coord[0] - pseudo_coord[0]
    dy = coord[1] - pseudo_coord[1]
    dz = coord[2] - pseudo_coord[2]
    distance = (dx**2 + dy**2 + dz**2) ** 0.5

# calculate LJ potential, may modify, in thesis, only **12 is considered, the ramaining 4 pi, - (radius/distance)**6 were omitted
    if distance < radius:
        lj_potential = (radius/distance)**12 
    else:
        lj_potential = 0.0

 # list_append, creating a list for individual results
    lj_potentials.append(lj_potential)

# sum LJ potentials
lj_sum = sum(lj_potentials)

# printing data finally
print("LJ potentials:")
for potential in lj_potentials:
    print(potential)  
