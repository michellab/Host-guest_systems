

import os
import re
import sys
import argparse
from Sire.Base import *
from datetime import datetime
# Make sure that the OPENMM_PLUGIN_DIR enviroment variable is set correctly.
os.environ["OPENMM_PLUGIN_DIR"] = getLibDir() + "/plugins"

from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Qt import *
from Sire.ID import *
from Sire.Config import *
from Sire.Analysis import *
from Sire.Tools.DCDFile import *
from Sire.Tools import Parameter, resolveParameters
import Sire.Stream
import time
import numpy as np
import mdtraj

perturbed_resnum = Parameter("perturbed residue number",1,"""The residue number of the molecule to morph.""")

morphfile = Parameter("morphfile", "../../input/MORPH.pert",
                      """Name of the morph file containing the perturbation to apply to the system.""")
cutoff_type = Parameter("cutoff type", "cutoffperiodic", """The cutoff method to use during the simulation.""")

cutoff_dist = Parameter("cutoff distance", 10 * angstrom,
                        """The cutoff distance to use for the non-bonded interactions.""")
rf_dielectric = Parameter("reaction field dielectric", 78.3,
                          """Dielectric constant to use if the reaction field cutoff method is used.""")

temperature = Parameter("temperature", 25 * celsius, """Simulation temperature""")

pressure = Parameter("pressure", 1 * atm, """Simulation pressure""")

lambda_val = Parameter("lambda_val", 0.0,
                       """Value of the lambda parameter at which to evaluate free energy gradients.""")

delta_lambda = Parameter("delta_lambda", 0.001,
                         """Value of the lambda interval used to evaluate free energy gradients by finite difference.""")

lambda_array = Parameter("lambda array",[] ,
                        """Array with all lambda values lambda_val needs to be part of the array. """)

shift_delta = Parameter("shift delta", 2.0,
                        """Value of the Lennard-Jones soft-core parameter.""")

coulomb_power = Parameter("coulomb power", 0,
                          """Value of the Coulombic soft-core parameter.""")

energy_frequency = Parameter("energy frequency", 1,
                             """The number of time steps between evaluation of free energy gradients.""")

# combining_rules = Parameter("combining rules", "geometric",
#                             """Combining rules to use for the non-bonded interactions.""")

timestep = Parameter("timestep", 2 * femtosecond, """Timestep for the dynamics simulation.""")

trajfile = Parameter("trajfile", "traj000000001.dcd",
                    """File name of the trajectory to process.""")
stepframe = Parameter("step_frame",1,"""The number of frames to step to between two succcessive evaluations.""")

def getDummies(molecule):
    print ("Selecting dummy groups")
    natoms = molecule.nAtoms()
    atoms = molecule.atoms()

    from_dummies = None
    to_dummies = None

    for x in range(0, natoms):
        atom = atoms[x]
        if atom.property("initial_ambertype") == "du":
            if from_dummies is None:
                from_dummies = molecule.selectAll(atom.index())
            else:
                from_dummies += molecule.selectAll(atom.index())
        elif atom.property("final_ambertype") == "du":
            if to_dummies is None:
                to_dummies = molecule.selectAll(atom.index())
            else:
                to_dummies += molecule.selectAll(atom.index())

    return to_dummies, from_dummies

def updateSystemfromTraj(system, frame_xyz, cell_lengths, cell_angles):
    traj_coordinates = frame_xyz[0]

    traj_box_x = float(cell_lengths[0][0])
    traj_box_y = float(cell_lengths[0][0])
    traj_box_z = float(cell_lengths[0][0])

    traj_natoms = len(traj_coordinates)

    # Sire does not support non rectangular boxes
    newmols_coords = {}

    traj_index = 0
    mol_index = 0

    molnums = system.molNums()
    molnums.sort()

    for molnum in molnums:
        mol = system.molecule(molnum)[0].molecule()
        molatoms = mol.atoms()
        molnatoms = mol.nAtoms()
        # Create an empty coord group using molecule so we get the correct layout
        newmol_coords = AtomCoords( mol.property("coordinates") )
        for x in range(0,molnatoms):
            tmparray = traj_coordinates[traj_index]
            atom_coord = Vector( tmparray[0].tolist() , tmparray[1].tolist() , tmparray[2].tolist() )
            atom = molatoms[x]
            cgatomidx = atom.cgAtomIdx()
            newmol_coords.set( cgatomidx, atom_coord)
            traj_index += 1
        newmols_coords[molnum] = newmol_coords
        mol_index += 1

    if traj_natoms != traj_index:
        print ("The number of atoms in the system is not equal to the number of atoms in the trajectory file ! Aborting.")
        sys.exit(-1)

    changedmols = MoleculeGroup("changedmols")
    mol_index = 0
    for molnum in molnums:
        mol = system.molecule(molnum)[0].molecule()
        newmol_coords = newmols_coords[molnum]
        mol = mol.edit().setProperty("coordinates", newmol_coords).commit()
        changedmols.add(mol)
    system.update(changedmols)

    # space = PeriodicBox(Vector( 79.7177, 79.7177, 79.7177 ) )
    space = PeriodicBox(Vector( traj_box_x, traj_box_y, traj_box_z ) )
    system.setProperty("space",space)

    return system


def createSystemFreeEnergy(molecules):
    r"""creates the system for free energy calculation
    Parameters
    ----------
    molecules : Sire.molecules
        Sire object that contains a lot of information about molecules
    Returns
    -------
    system : Sire.system
    """
    print ("Create the System...")

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber)[0].molecule()
        moleculeList.append(molecule)

    # Scan input to find a molecule with passed residue number 
    # The residue name of the first residue in this molecule is
    # used to name the solute. This is used later to match
    # templates in the flex/pert files.

    solute = None
    for molecule in moleculeList:
        if ( molecule.residue(ResIdx(0)).number() == ResNum(perturbed_resnum.val) ):
            solute = molecule
            moleculeList.remove(molecule)
            break

    if solute is None:
        print ("FATAL ! Could not find a solute to perturb with residue number %s in the input ! Check the value of your cfg keyword 'perturbed residue number'" % perturbed_resnum.val)
        sys.exit(-1)

    #solute = moleculeList[0]
    host_1 =  moleculeList[0]
    host_2 =  moleculeList[1]
    host_3 =  moleculeList[2]
    host_4 =  moleculeList[3]
    pd1 = moleculeList[4]
    pd2 = moleculeList[5]
    ion1 = moleculeList[6]
    ion2 = moleculeList[7]
    ion3 = moleculeList[8]
    ion4 = moleculeList[9]

    solv = []
    for i in range(10, len(moleculeList)):
        solv.append(moleculeList[i])

    lig_name = solute.residue(ResIdx(0)).name().value()

    solute = solute.edit().rename(lig_name).commit()

    perturbations_lib = PerturbationsLibrary(morphfile.val)
    solute = perturbations_lib.applyTemplate(solute)

    perturbations = solute.property("perturbations")

    lam = Symbol("lambda")

    initial = Perturbation.symbols().initial()
    final = Perturbation.symbols().final()

    solute = solute.edit().setProperty("perturbations",
                                       perturbations.recreate((1 - lam) * initial + lam * final)).commit()

    # We put atoms in three groups depending on what happens in the perturbation
    # non dummy to non dummy --> the hard group, uses a normal intermolecular FF
    # non dummy to dummy --> the todummy group, uses SoftFF with alpha = Lambda
    # dummy to non dummy --> the fromdummy group, uses SoftFF with alpha = 1 - Lambda
    # We start assuming all atoms are hard atoms. Then we call getDummies to find which atoms
    # start/end as dummies and update the hard, todummy and fromdummy groups accordingly

    solute_grp_ref = MoleculeGroup("solute_ref", solute)

    # In [27]: type(solute_grp_ref.molecules())                                                                                
    # Out[27]: Sire.Mol._Mol.Molecules

    # In [28]: type(solute_grp_ref.molecules().first())                                                                        
    # Out[28]: Sire.Mol._Mol.Molecule

    # these are empty. no molecules in the set 
    solute_grp_ref_hard = MoleculeGroup("solute_ref_hard")
    solute_grp_ref_todummy = MoleculeGroup("solute_ref_todummy")
    solute_grp_ref_fromdummy = MoleculeGroup("solute_ref_fromdummy")


    solute_ref_hard = solute.selectAllAtoms()
    solute_ref_todummy = solute_ref_hard.invert()
    solute_ref_fromdummy = solute_ref_hard.invert()

    # create a new group with the host, ions and solvent
    host1_hard = host_1.selectAllAtoms()
    host2_hard = host_2.selectAllAtoms()
    host3_hard = host_3.selectAllAtoms()
    host4_hard = host_4.selectAllAtoms()
    pd1_hard = pd1.selectAllAtoms()
    pd2_hard = pd2.selectAllAtoms()
    ion1_hard = ion1.selectAllAtoms()
    ion2_hard = ion2.selectAllAtoms()
    ion3_hard = ion3.selectAllAtoms()
    ion4_hard = ion4.selectAllAtoms()


    # solv_hard = solv.selectAllAtoms() 

    host1_grp_hard = MoleculeGroup("host1_hard")
    # host1_hard = host_1.selectAllAtoms()
    host1_grp_hard.add(host1_hard)

    host_ions_grp_hard = MoleculeGroup("host_ions_hard")
    host_ions_grp_hard.add(host1_hard)
    host_ions_grp_hard.add(host2_hard)
    host_ions_grp_hard.add(host3_hard)
    host_ions_grp_hard.add(host4_hard)
    host_ions_grp_hard.add(pd1_hard)
    host_ions_grp_hard.add(pd2_hard)
    host_ions_grp_hard.add(ion1_hard)
    host_ions_grp_hard.add(ion2_hard)
    host_ions_grp_hard.add(ion3_hard)
    host_ions_grp_hard.add(ion4_hard)

    host_ions_solv_grp_hard = MoleculeGroup("host_ions_solvent_hard")
    host_ions_solv_grp_hard.add(host1_hard)
    host_ions_solv_grp_hard.add(host2_hard)
    host_ions_solv_grp_hard.add(host3_hard)
    host_ions_solv_grp_hard.add(host4_hard)
    host_ions_solv_grp_hard.add(pd1_hard)
    host_ions_solv_grp_hard.add(pd2_hard)
    host_ions_solv_grp_hard.add(ion1_hard)
    host_ions_solv_grp_hard.add(ion2_hard)
    host_ions_solv_grp_hard.add(ion3_hard)
    host_ions_solv_grp_hard.add(ion4_hard)
    for i in range (0, len(solv)):
        host_ions_solv_grp_hard.add(solv[i]) 


    to_dummies, from_dummies = getDummies(solute)

    if to_dummies is not None:
        ndummies = to_dummies.count()
        dummies = to_dummies.atoms()

        for x in range(0, ndummies):
            dummy_index = dummies[x].index()
            solute_ref_hard = solute_ref_hard.subtract(solute.select(dummy_index))
            solute_ref_todummy = solute_ref_todummy.add(solute.select(dummy_index))

    if from_dummies is not None:
        ndummies = from_dummies.count()
        dummies = from_dummies.atoms()

        for x in range(0, ndummies):
            dummy_index = dummies[x].index()
            solute_ref_hard = solute_ref_hard.subtract(solute.select(dummy_index))
            solute_ref_fromdummy = solute_ref_fromdummy.add(solute.select(dummy_index))

    solute_grp_ref_hard.add(solute_ref_hard)
    solute_grp_ref_todummy.add(solute_ref_todummy)
    solute_grp_ref_fromdummy.add(solute_ref_fromdummy)

    solutes = MoleculeGroup("solutes")
    solutes.add(solute)

    molecules = MoleculeGroup("molecules")
    molecules.add(solute)

    solvent = MoleculeGroup("solvent")

    #for molecule in moleculeList[1:]:
    for molecule in moleculeList:
        molecules.add(molecule)

    # In [9]: molecules                                                                                                        
    # Out[9]: MoleculeGroup(molecules, 13)

    # In [10]: molecules.molecules()                                                                                           
    # Out[10]: Molecules{ nMolecules() == 4436, nViews() == 4436 }

    for molecule in solv:
        solvent.add(molecule)

    # In [7]: solvent                                                                                                          
    # Out[7]: MoleculeGroup(solvent, 14)

    # In [8]: solvent.molecules()                                                                                              
    # Out[8]: Molecules{ nMolecules() == 4425, nViews() == 4425 }
    non_guest = MoleculeGroup("non_guest")

    #for molecule in moleculeList[1:]:
    for molecule in moleculeList:
        # molecules.add(molecule)
        non_guest.add(molecule)



    all = MoleculeGroup("all")

    all.add(molecules)
    all.add(solvent)

    # "non-guest" is the old "solvent" that contains everything but the perturbed molecule
    # Will use it at the calculation of the total energy
    all.add(non_guest)

    all.add(solutes)
    all.add(solute_grp_ref)
    all.add(solute_grp_ref_hard)
    all.add(solute_grp_ref_todummy)
    all.add(solute_grp_ref_fromdummy)
    all.add(host_ions_grp_hard)
    all.add(host_ions_solv_grp_hard)

    # Add these groups to the System
    system = System()

    system.add(solutes)
    system.add(solute_grp_ref)
    system.add(solute_grp_ref_hard)
    system.add(solute_grp_ref_todummy)
    system.add(solute_grp_ref_fromdummy)

    system.add(molecules)

    system.add(solvent)
    system.add(non_guest)

    system.add(host_ions_grp_hard)
    system.add(host_ions_solv_grp_hard)

    system.add(all)

    return system


def setupForceFieldsFreeEnergy(system, space):
    r"""sets up the force field for the free energy calculation
    Parameters
    ----------
    system : Sire.system
    space : Sire.space
    Returns
    -------
    system : Sire.system
    """

    print ("Creating force fields... ")

    solutes = system[MGName("solutes")]

    solute = system[MGName("solute_ref")]
    solute_hard = system[MGName("solute_ref_hard")]
    solute_todummy = system[MGName("solute_ref_todummy")]
    solute_fromdummy = system[MGName("solute_ref_fromdummy")]

    host_ions_solv_grp_hard = system[MGName("host_ions_solvent_hard")]
    host_ions_grp_hard = system[MGName("host_ions_hard")]

    solvent = system[MGName("solvent")]
    non_guest = system[MGName("non_guest")]


    all = system[MGName("all")]

    # ''solvent'' is just the dichloromethane molecules
    solvent_intraff = InternalFF("solvent_intraff")
    solvent_intraff.add(solvent)

    non_guest_intraff = InternalFF("non_guest_intraff")
    non_guest_intraff.add(non_guest)

    # Guest bond, angle, dihedral energy
    solute_intraff = InternalFF("solute_intraff")
    solute_intraff.add(solute)

    # host and ions  bond, angle, dihedral energy
    host_ions_intraff = InternalFF("host_ions_intraff")
    host_ions_intraff.add(host_ions_grp_hard)

    # host, ions and solvent bond, angle, dihedral energy
    host_ions_solvent_intraff = InternalFF("host_ions_solvent_intraff")
    host_ions_solvent_intraff.add(host_ions_solv_grp_hard)

    # Solvent-solvent coulomb/LJ (CLJ) energy
    solventff = InterCLJFF("solvent:solvent")
    if (cutoff_type.val != "nocutoff"):
        solventff.setUseReactionField(True)
        solventff.setReactionFieldDielectric(rf_dielectric.val)
    solventff.add(solvent)

    # everything but guest coulomb/LJ (CLJ) energy
    non_guestff = InterCLJFF("non_guest:non_guest")
    if (cutoff_type.val != "nocutoff"):
        non_guestff.setUseReactionField(True)
        non_guestff.setReactionFieldDielectric(rf_dielectric.val)
    non_guestff.add(non_guest)

    #Solvent intramolecular CLJ energy
    solvent_intraclj = IntraCLJFF("solvent_intraclj")
    if (cutoff_type.val != "nocutoff"):
        solvent_intraclj.setUseReactionField(True)
        solvent_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solvent_intraclj.add(solvent)

    #everything but guest intramolecular CLJ energy
    non_guest_intraclj = IntraCLJFF("non_guest_intraclj")
    if (cutoff_type.val != "nocutoff"):
        non_guest_intraclj.setUseReactionField(True)
        non_guest_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    non_guest_intraclj.add(non_guest)

    # host, ions, solvent intramolecular CLJ energy
    host_ions_solvent_hard_intraclj = IntraCLJFF("host_ions_solvent_hard_intraclj")
    if (cutoff_type.val != "nocutoff"):
        host_ions_solvent_hard_intraclj.setUseReactionField(True)
        host_ions_solvent_hard_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    host_ions_solvent_hard_intraclj.add(host_ions_solv_grp_hard)

    # host, ions intramolecular CLJ energy
    host_ions_hard_intraclj = IntraCLJFF("host_ions_hard_intraclj")
    if (cutoff_type.val != "nocutoff"):
        host_ions_hard_intraclj.setUseReactionField(True)
        host_ions_hard_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    host_ions_hard_intraclj.add(host_ions_grp_hard)

    # Solute intramolecular CLJ energy
    solute_hard_intraclj = IntraCLJFF("solute_hard_intraclj")
    if (cutoff_type.val != "nocutoff"):
        solute_hard_intraclj.setUseReactionField(True)
        solute_hard_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_intraclj.add(solute_hard)

    solute_todummy_intraclj = IntraSoftCLJFF("solute_todummy_intraclj")
    solute_todummy_intraclj.setShiftDelta(shift_delta.val)
    solute_todummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_intraclj.setUseReactionField(True)
        solute_todummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_intraclj.add(solute_todummy)

    solute_fromdummy_intraclj = IntraSoftCLJFF("solute_fromdummy_intraclj")
    solute_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    solute_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_fromdummy_intraclj.setUseReactionField(True)
        solute_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_fromdummy_intraclj.add(solute_fromdummy)

    solute_hard_todummy_intraclj = IntraGroupSoftCLJFF("solute_hard:todummy_intraclj")
    solute_hard_todummy_intraclj.setShiftDelta(shift_delta.val)
    solute_hard_todummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_hard_todummy_intraclj.setUseReactionField(True)
        solute_hard_todummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_todummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_todummy_intraclj.add(solute_todummy, MGIdx(1))

    solute_hard_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_hard:fromdummy_intraclj")
    solute_hard_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    solute_hard_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_hard_fromdummy_intraclj.setUseReactionField(True)
        solute_hard_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_fromdummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))

    solute_todummy_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_todummy:fromdummy_intraclj")
    solute_todummy_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    solute_todummy_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_fromdummy_intraclj.setUseReactionField(True)
        solute_todummy_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_fromdummy_intraclj.add(solute_todummy, MGIdx(0))
    solute_todummy_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))

    #Solute-guest, ions, solvent  CLJ energy
    solute_guest_ions_solvent_hard_solventff = InterGroupCLJFF("solute_hard:guest-ions-solvent")
    if (cutoff_type.val != "nocutoff"):
        solute_guest_ions_solvent_hard_solventff.setUseReactionField(True)
        solute_guest_ions_solvent_hard_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_guest_ions_solvent_hard_solventff.add(solute_hard, MGIdx(0))
    solute_guest_ions_solvent_hard_solventff.add(host_ions_solv_grp_hard , MGIdx(1))

    solute_guest_ions_solvent_todummy_solventff = InterGroupSoftCLJFF("solute_todummy:guest-ions-solvent")
    if (cutoff_type.val != "nocutoff"):
        solute_guest_ions_solvent_todummy_solventff.setUseReactionField(True)
        solute_guest_ions_solvent_todummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_guest_ions_solvent_todummy_solventff.add(solute_todummy, MGIdx(0))
    solute_guest_ions_solvent_todummy_solventff.add(host_ions_solv_grp_hard, MGIdx(1))

    solute_guest_ions_solvent_fromdummy_solventff = InterGroupSoftCLJFF("solute_fromdummy:guest-ions-solvent")
    if (cutoff_type.val != "nocutoff"):
        solute_guest_ions_solvent_fromdummy_solventff.setUseReactionField(True)
        solute_guest_ions_solvent_fromdummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_guest_ions_solvent_fromdummy_solventff.add(solute_fromdummy, MGIdx(0))
    solute_guest_ions_solvent_fromdummy_solventff.add(host_ions_solv_grp_hard, MGIdx(1))

    #Solute-solvent CLJ energy
    solute_hard_solventff = InterGroupCLJFF("solute_hard:solvent")
    if (cutoff_type.val != "nocutoff"):
        solute_hard_solventff.setUseReactionField(True)
        solute_hard_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_solventff.add(solute_hard, MGIdx(0))
    solute_hard_solventff.add(solvent, MGIdx(1))

    solute_todummy_solventff = InterGroupSoftCLJFF("solute_todummy:solvent")
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_solventff.setUseReactionField(True)
        solute_todummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_solventff.add(solute_todummy, MGIdx(0))
    solute_todummy_solventff.add(solvent, MGIdx(1))

    solute_fromdummy_solventff = InterGroupSoftCLJFF("solute_fromdummy:solvent")
    if (cutoff_type.val != "nocutoff"):
        solute_fromdummy_solventff.setUseReactionField(True)
        solute_fromdummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_fromdummy_solventff.add(solute_fromdummy, MGIdx(0))
    solute_fromdummy_solventff.add(solvent, MGIdx(1))

    #non_guest-non_guest CLJ energy
    solute_hard_non_guestff = InterGroupCLJFF("solute_hard:non_guest")
    if (cutoff_type.val != "nocutoff"):
        solute_hard_non_guestff.setUseReactionField(True)
        solute_hard_non_guestff.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_non_guestff.add(solute_hard, MGIdx(0))
    solute_hard_non_guestff.add(non_guest, MGIdx(1))

    solute_todummy_non_guestff = InterGroupSoftCLJFF("solute_todummy:non_guest")
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_non_guestff.setUseReactionField(True)
        solute_todummy_non_guestff.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_non_guestff.add(solute_todummy, MGIdx(0))
    solute_todummy_non_guestff.add(non_guest, MGIdx(1))

    solute_fromdummy_non_guestff = InterGroupSoftCLJFF("solute_fromdummy:non_guest")
    if (cutoff_type.val != "nocutoff"):
        solute_fromdummy_non_guestff.setUseReactionField(True)
        solute_fromdummy_non_guestff.setReactionFieldDielectric(rf_dielectric.val)
    solute_fromdummy_non_guestff.add(solute_fromdummy, MGIdx(0))
    solute_fromdummy_non_guestff.add(non_guest, MGIdx(1))


    # TOTAL
    forcefields = [solute_intraff,
                   solute_hard_intraclj, solute_todummy_intraclj, solute_fromdummy_intraclj,
                   host_ions_solvent_hard_intraclj, solute_hard_todummy_intraclj,
                    solute_hard_fromdummy_intraclj, solute_todummy_fromdummy_intraclj,
                   solute_guest_ions_solvent_hard_solventff, solute_guest_ions_solvent_fromdummy_solventff,
                   solute_guest_ions_solvent_todummy_solventff, host_ions_hard_intraclj,
                   solvent_intraff,
                   solventff, solvent_intraclj,
                   solute_hard_solventff, solute_todummy_solventff, solute_fromdummy_solventff,
                   non_guest_intraff,
                   non_guestff, non_guest_intraclj,
                   solute_hard_non_guestff, solute_todummy_non_guestff, solute_fromdummy_non_guestff]


    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)

    if (cutoff_type.val != "nocutoff"):
        system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff_dist.val))
    else:
        system.setProperty("switchingFunction", NoCutoff())

    # system.setProperty("combiningRules", VariantProperty(combining_rules.val))
    system.setProperty("coulombPower", VariantProperty(coulomb_power.val))
    system.setProperty("shiftDelta", VariantProperty(shift_delta.val))

    # TOTAL
    total_nrg = solute_intraff.components().total() + solute_hard_intraclj.components().total() + \
                solute_todummy_intraclj.components().total(0) + solute_fromdummy_intraclj.components().total(0) + \
                solute_hard_todummy_intraclj.components().total(
                    0) + solute_hard_fromdummy_intraclj.components().total(0) + \
                solute_todummy_fromdummy_intraclj.components().total(0) + \
                non_guest_intraff.components().total() + non_guestff.components().total() + \
                non_guest_intraclj.components().total() + \
                solute_hard_non_guestff.components().total() + \
                solute_todummy_non_guestff.components().total(0) + \
                solute_fromdummy_non_guestff.components().total(0)

    e_total = system.totalComponent()

    lam = Symbol("lambda")

    system.setComponent(e_total, total_nrg)

    system.setConstant(lam, 0.0)

    system.add(PerturbationConstraint(solutes))

    # NON BONDED Alpha constraints for the soft force fields

    system.add(PropertyConstraint("alpha0", FFName("solute_todummy_intraclj"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_fromdummy_intraclj"), 1 - lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_hard:todummy_intraclj"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_hard:fromdummy_intraclj"), 1 - lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_todummy:fromdummy_intraclj"), Max(lam, 1 - lam)))
    system.add(PropertyConstraint("alpha0", FFName("solute_todummy:solvent"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_fromdummy:solvent"), 1 - lam))

    system.setComponent(lam, lambda_val.val)

    # printEnergies( system.componentValues() )

    return system


crdfile = "../../input/SYSTEM.crd"
topfile = "../../input/SYSTEM.top"

amber = Amber()
(molecules, space) = amber.readCrdTop(crdfile, topfile)
system = createSystemFreeEnergy(molecules)
system = setupForceFieldsFreeEnergy(system, space)

# Now loop over snapshots in dcd, discard the first 20%
start_frame = 0
end_frame = 10
step_frame = stepframe.val

#mdtraj_top = mdtraj.load_prmtop(topfile.val)
mdtraj_trajfile = mdtraj.open(trajfile.val,'r')
nframes = len(mdtraj_trajfile)
if end_frame > (nframes - 1):
    end_frame = nframes - 1
mdtraj_trajfile.seek(start_frame)
current_frame = start_frame

sys_nrgs = []

while (current_frame <= end_frame):
    print ("Processing frame %s " % current_frame)
    frames_xyz, cell_lengths, cell_angles = mdtraj_trajfile.read(n_frames=1)
    system = updateSystemfromTraj(system, frames_xyz, cell_lengths, cell_angles)
    sys_nrg = system.energies()
    sys_nrgs.append(sys_nrg)
    current_frame += step_frame
    mdtraj_trajfile.seek(current_frame)
#print (sys_nrgs)

# open a new file to write
outF = open("h-g_energies.txt", "w")
md_str = str(sys_nrgs)
outF.write(md_str)
outF.close()

#mdtraj_top = mdtraj.load_prmtop(topfile.val)
mdtraj_trajfile = mdtraj.open(trajfile.val,'r')
nframes = len(mdtraj_trajfile)
if end_frame > (nframes - 1):
    end_frame = nframes - 1
mdtraj_trajfile.seek(start_frame)
current_frame = start_frame

sum_sys_nrg = 0
for i in range(1, end_frame):
    sum_sys_nrg = sum_sys_nrg + sys_nrgs[i].value()
avegare_sys_nrg = sum_sys_nrg/(end_frame - start_frame)
print(avegare_sys_nrg, "kcal/mol")
