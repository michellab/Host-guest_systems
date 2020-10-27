# Host-guest systems

Notes: The final amber files for the combined system can be found at amber_files/combined_system, after optimising the protocol to get the correct geometries for the ligands (see below for details).

Protocol: 
* host ligands: Parameters generated from LigParGen, amber files from qube_to_amber.py
* guest: Parameters generated from LigParGen, amber files from qube_to_amber.py
* Pd(II): Parameters generated from Macromodel, amber files from BSS conversion from gromacs files
* barf: Parameters generated from Macromodel, amber files from BSS conversion from gromacs files
* dichloromethane:  Parameters generated from Macromodel, amber files from BSS conversion from gromacs files

The coordinates of the molecules were extracted from a solvated system using gromacs. For the LigParGen molecules (guest and host ligands), the pdb files had different atom typed, so the LigParGen pdb files were aligned with the gromacs files and the correct coordinated were saved (initial_conf_pdb_files). Then all files were converted to amber-type ones (amber_files/initial_structures) and combined to one system (amber_files/combined_system) using BSS. Then MORPH.pert files were generated for the discharge and vanish steps using morph_step*.py.  

## Conversions to amber files: 
a) Using BSS 
```
~/biosimm2019.app/bin/ipython
import BioSimSpace as BSS
BSS.setVerbose(True)
files = ["conf.gro","topol.top"]
system = BSS.IO.readMolecules(files, property_map={"GROMACS_PATH" : "/home/sofia/Desktop/Sofia/cages/LigParGen_topol/top"})
BSS.IO.saveMolecules("OUTPUT", system, ["rst7", "prm7"])
quit()
```

b) Using qube_to_amber.py with XML files from LigParGen.

## Comparison of MD trajectories of the guest G4 using GROMACS and SOMD:
- GROMACS: Structure of the ligand is correct. The molecule doesn't bend. 
- SOMD with files from BSS conversion: Ligands don't maintain the correct geometry. Problem with dihedrals as aromatic moieties bend. 
- SOMD with files from LigParGen and qube_to_amber.py conversion: Ligands have the correct geometry. No obvious issues.

## Energy comparison with Sire unit test:
https://github.com/michellab/SireUnitTests/blob/devel/unittests/SireIO/test_groamb.py

Dihedrals are not converted correctly. Also see issue here: https://github.com/michellab/BioSimSpace/issues/120

```
E_{cljff}^{coulomb_{default}}:  -4285.729874817867  versus  -4285.729874817867
E_{intraclj}^{CLJ_{default}}:  -625.2446171181026  versus  -625.2446171181026
E_{intraclj}^{LJ_{default}}:  -24.58456319891593  versus  -24.58456319891593
E_{intraclj}^{coulomb_{default}}:  -600.6600539191867  versus  -600.6600539191867
E_{intraff}^{1-4[LJ]}:  257.8395504590764  versus  257.8395512628246
E_{intraff}^{1-4[coulomb]}:  297.32472276315536  versus  297.32472276315536
E_{intraff}^{1-4}:  555.1642732222317  versus  555.1642740259799
E_{intraff}^{angle}:  7559.450773754545  versus  7559.450565986063
E_{intraff}^{bend-bend}:  0.0  versus  0.0
E_{intraff}^{bond}:  88.12330615066351  versus  88.12330615068086
E_{intraff}^{dihedral}:  78.51547305041994  versus  8.953272101907439
E_{intraff}^{improper}:  0.0  versus  0.0
E_{intraff}^{internal}:  8281.25382617786  versus  8211.691418264632
E_{intraff}^{stretch-bend-torsion}:  0.0  versus  0.0
E_{intraff}^{stretch-bend}:  0.0  versus  0.0
E_{intraff}^{stretch-stretch}:  0.0  versus  0.0
E_{intraff}^{urey_bradley}:  0.0  versus  0.0
```
## Parameters:

To maintain the correct geometry of the **host ligands** and the **guest**, the LigParGen xml files were used and converted to amber files with qube_to_amber.py
The counter ion **BARF-** and **Pd(II)** could not be parameterised with LigParGen, so the Macromodel OPLS-AA parameters were used instead. 
**Dichloromethane** was used as solvent with parameters from Macromodel that were optimised using non bonded parameters from Scheurer _et al._ https://pubs.acs.org/doi/pdf/10.1021/acs.jpcb.7b06709

