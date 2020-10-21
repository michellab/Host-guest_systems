The combination of host/guest molecules with the Pd(II) ions to form the capsule and the bound system are done with VMD, using TopoTools at VMD's TK console. 
Firstly the molecules are moved/rotated and saved to get the correct positions. Then we append a list with all the molecules and we use the VMD plugin "TopoTools" to write a pdb with the combined system.

The protocol at "protocol.txt" uses GROMACS to solvate the molecular system with dichloromethane and then it minimises and equilibrates it.
- Minimisation: 
```
emtol  = 100.0         ; Stop minimization when the maximum force < 100.0 kJ/mol
nsteps = 50000         ; Maximum number of (minimization) steps to perform
```
                
- NVT equilibration: 
```
nsteps = 50000    ; 2 * 50000 = 100 ps
dt     = 0.002    ; 2 fs
```

- NPT equilibration: 
```
nsteps = 50000    ; 2 * 50000 = 100 ps
dt     = 0.002    ; 2 fs
```
                     
