This protocol uses GROMACS to solvate the molecular system with dichloromethane and then it minimises and equilibrates it.
- Minimisation: emtol  = 100.0         ; Stop minimization when the maximum force < 100.0 kJ/mol
                nsteps = 50000         ; Maximum number of (minimization) steps to perform
                
- NVT equilibration: nsteps = 50000    ; 2 * 50000 = 100 ps
                     dt     = 0.002    ; 2 fs

- NPT equilibration: nsteps = 50000    ; 2 * 50000 = 100 ps
                     dt     = 0.002    ; 2 fs
                     
