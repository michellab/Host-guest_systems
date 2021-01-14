Description of files: 

- **ener-over-snapshots.py** 

Calculates the energy components over every snapshot and writes them in a file. 
Parameters to change are the first and the last frame of the trajectory. Ususally discard the first 20% of the frames

- **host-guest_ener_over_snapshots.py**

Same as ener-over-snapshots.py, but with  new groups added for the calculation of the energetic contributions: 
host-ions, host-ions-solvent and non_guest
Other modifications: 
i) solvent refers to all the solvent molecules (dichloromethane), and not to all the non-perturbed molecules as ener-over-snapshots.py does
ii) to get the same total potential energy the group "non_guest" is added, that includes all the non-perturbed molecules - same as solvent of ener-over-snapshots.py
Indeed the two scripts give the same total energy: 
```
output from energies.txt:  E_{total} == -17184.8 
output from h-g_energies.txt:  E_{total} == -17184.8 
```

- **scale-charges.py**

Exports a new topology file with modified charges. New charges = old_charges/ scale_constant 
Usage: `python scale_charges.py  scale_constant num_of_atoms`

- **calc_average_ener.py**

Returns averages of the energetic components 
