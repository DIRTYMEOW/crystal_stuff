This is a brief and perhaps oversimplified estimation of vdw interaction contributed from neighboring atoms to a "Ghost atom".
The position of ghost atom, and xyz file has to be prepared first.
Define a coordination of ghost atom, modify the .py file.
.xyz file name is here~~~~~~~~~~~~~~~ for other file, u need to re-write coding here.
The 12-6 equation, u can tune whatever value u want, or even vdw radii.
Caution:
This calculation overlooks the concept of "molecules" "orbital" "dipole" "charge" "bonding", just treat everthing as single atom, which acts simply as a lego.

The original purpose was to estimate the possibility of photoisomerization of single molecules within single crystal,
By generating a n*n*n crystal from cif file into .xyz file, I used it to estimate the repulsion contributed at a single point of molecule,
and see if those having lower vdw repulsion are photoisomerizable, and vice versa.

Feel free to contact if u need further help.
d05223110@ntu.edu.tw
~~~~~~~~~~~~~~
First, "mkdir big_xyz pedaling_ghost single_crystal"
Put all ur *.cif in ./single_crystal
run "python step1_grep_cif_to_xyz.py", in this file, u can modify how many n*n*n u want.
Chk the file generated in the ./big_xyz, see if the molecules are many enough, normally 10*10*10 will do.
run "python step2_replace_central_to_ghost.py", this code identify the most central molecule, and replace it into a ghost atom (Barium in the code), in the original file, I identify the central molecule by recognizing C=N (C-N bond shorter than 1.3 ang.). This part could be tricky, if u need help, can contact me or use ChatGPT to modify. After successful running, chk ./pedaling_ghost, to see if the central molecule is replaced into a barium atom.
run "python step3_vdw_dirty_energy", where u can finetune the radius of ghost_atom, which is Ba here.
In step3-2_dirty_space.py, the code is still scratchy, but it estimate void around the Ba. U can finetune the fixed_radius = x, so the range of spherical of Ba can be assigned. The occupied space is based on vdw surface generated from the nearest molecules.
