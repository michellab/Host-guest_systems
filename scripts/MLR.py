
import os
import pandas as pd
import numpy as np 
import statsmodels.api as sm


def hard_sol_C(): 
    hard_solv_C = []
    count = 0
    iname=("h-g_energies.txt")
    ener_list = []
    with open(iname) as infile:
        lines = infile.readlines()
        for line in lines:
            for i in range (0, len(line.split()) ):
                if line.split()[i] == "E_{solute_hard:solvent}^{coulomb}":
                    hard_solv_C.append(line.split()[i+2].strip(","))
                    entry = float(line.split()[i+2].strip(","))
                    count = count + entry
    # print ("hard_solv_C = ",hard_solv_C)
    avg = count/len(hard_solv_C)
    print("Average E_{solute_hard:solvent}^{coulomb} = ", avg)
    return(avg)


def hard_sol_LJ(): 
    hard_solv_LJ = []
    count = 0
    iname=("h-g_energies.txt")
    ener_list = []
    with open(iname) as infile:
        lines = infile.readlines()
        for line in lines:
            for i in range (0, len(line.split()) ):
                if line.split()[i] == "E_{solute_hard:solvent}^{LJ}":
                    hard_solv_LJ.append(line.split()[i+2].strip(","))
                    entry = float(line.split()[i+2].strip(","))
                    count = count + entry
    # print ("hard_solv_LJ = ",hard_solv_LJ)
    avg = count/len(hard_solv_LJ)
    print("Average E_{solute_hard:solvent}^{LJ} = ", avg)
    return(avg)



def hard_intraC(): 
    hard_intra_C = []
    count = 0
    iname=("h-g_energies.txt")
    ener_list = []
    with open(iname) as infile:
        lines = infile.readlines()
        for line in lines:
            for i in range (0, len(line.split()) ):
                if line.split()[i] == "E_{solute_hard_intraclj}^{coulomb}":
                    hard_intra_C.append(line.split()[i+2].strip(","))
                    entry = float(line.split()[i+2].strip(","))
                    count = count + entry
    # print ("hard_intra_C = ",hard_intra_C)
    avg = count/len(hard_intra_C)
    print("Average E_{solute_hard_intraclj}^{coulomb} = ", avg)
    return(avg)


def hard_intraLJ(): 
    hard_intra_LJ = []
    count = 0
    iname=("h-g_energies.txt")
    ener_list = []
    with open(iname) as infile:
        lines = infile.readlines()
        for line in lines:
            for i in range (0, len(line.split()) ):
                if line.split()[i] == "E_{solute_hard_intraclj}^{LJ}":
                    hard_intra_LJ.append(line.split()[i+2].strip(","))
                    entry = float(line.split()[i+2].strip(","))
                    count = count + entry
    # print ("hard_intra_LJ = ",hard_intra_LJ)
    avg = count/len(hard_intra_LJ)
    print("Average E_{solute_hard_intraclj}^{LJ}J = ", avg)
    return(avg)


def hard_intraff14C(): 
    hard_intraff14_C = []
    count = 0
    iname=("h-g_energies.txt")
    ener_list = []
    with open(iname) as infile:
        lines = infile.readlines()
        for line in lines:
            for i in range (0, len(line.split()) ):
                if line.split()[i] == "E_{solute_intraff}^{1-4[coulomb]}":
                    hard_intraff14_C.append(line.split()[i+2].strip(","))
                    entry = float(line.split()[i+2].strip(","))
                    count = count + entry
    # print ("hard_intraff14_C = ",hard_intraff14_C)
    avg = count/len(hard_intraff14_C)
    print("Average E_{solute_intraff}^{1-4[coulomb]} = ", avg)
    return(avg)


def hard_intraff14LJ(): 
    hard_intraff14_LJ = []
    count = 0
    iname=("h-g_energies.txt")
    ener_list = []
    with open(iname) as infile:
        lines = infile.readlines()
        for line in lines:
            for i in range (0, len(line.split()) ):
                if line.split()[i] == "E_{solute_intraff}^{1-4[LJ]}":
                    hard_intraff14_LJ.append(line.split()[i+2].strip(","))
                    entry = float(line.split()[i+2].strip(","))
                    count = count + entry
    # print ("hard_intraff14_LJ = ",hard_intraff14_LJ)
    avg = count/len(hard_intraff14_LJ)
    print("E_{solute_intraff}^{1-4[LJ]} = ", avg)
    return(avg)



def sol_sol_LJ(): 
    solv_solv_LJ = []
    count = 0
    iname=("h-g_energies.txt")
    ener_list = []
    with open(iname) as infile:
        lines = infile.readlines()
        for line in lines:
            for i in range (0, len(line.split()) ):
                if line.split()[i] == "E_{solvent:solvent}^{LJ}":
                    solv_solv_LJ.append(line.split()[i+2].strip(","))
                    entry = float(line.split()[i+2].strip(","))
                    count = count + entry
    # print ("solv_solv_LJ = ",solv_solv_LJ)
    avg = count/len(solv_solv_LJ)
    print("Average E_{solvent:solvent}^{LJ} = ", avg)
    return(avg)

       
def sol_sol_C(): 
    solv_solv_C = []
    count = 0
    iname=("h-g_energies.txt")
    ener_list = []
    with open(iname) as infile:
        lines = infile.readlines()
        for line in lines:
            for i in range (0, len(line.split()) ):
                if line.split()[i] == "E_{solvent:solvent}^{coulomb}":
                    solv_solv_C.append(line.split()[i+2].strip(","))
                    entry = float(line.split()[i+2].strip(","))
                    count = count + entry
    # print ("solv_solv_C = ",solv_solv_C)
    avg = count/len(solv_solv_C)
    print("Average E_{solvent:solvent}^{coulomb} = ", avg)
    return(avg)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cage_in_solvent_elec = -68.6378272
cage_in_solvent_LJ = -32.384


list_E_free_C = []
list_E_free_LJ = []

print("---------------- Cage G1 ----------------")

os.chdir('/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G1/free/discharge/output/lambda-0.00')

hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul

hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs
sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()


E_free_C = total_host_guest_Coul + sol_sol_Coul - sol_sol_Coul
E_free_LJ = total_host_guest_LJs + sol_sol_LJs - sol_sol_LJs

list_E_free_C.append(E_free_C)
list_E_free_LJ.append(E_free_LJ)

print("list_E_free_C= ", list_E_free_C)
print("list_E_free_LJ= ", list_E_free_LJ)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

os.chdir('/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G1/free/discharge/output/lambda-0.00')
print("working direcory =", os.getcwd())
list_E_bound_C = []
list_E_bound_LJ = []


hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul
hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs

sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()

E_bound_C = total_host_guest_Coul  - cage_in_solvent_elec
E_bound_LJ = total_host_guest_LJs  - cage_in_solvent_LJ

list_E_bound_C.append(E_bound_C)
list_E_bound_LJ.append(E_bound_LJ)

print("list_E_bound_C= ", list_E_bound_C)
print("list_E_bound_LJ= ", list_E_bound_LJ)

os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G2/original/bound/discharge/output/lambda-0.00')

print("---------------- Cage G2 ----------------")


hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul

hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs
sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()


E_free_C = total_host_guest_Coul + sol_sol_Coul - sol_sol_Coul
E_free_LJ = total_host_guest_LJs + sol_sol_LJs - sol_sol_LJs

list_E_free_C.append(E_free_C)
list_E_free_LJ.append(E_free_LJ)

print("list_E_free_C= ", list_E_free_C)
print("list_E_free_LJ= ", list_E_free_LJ)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G2/original/free/discharge/output/lambda-0.00')
print("working direcory =", os.getcwd())

hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul
hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs

sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()

E_bound_C = total_host_guest_Coul  - cage_in_solvent_elec
E_bound_LJ = total_host_guest_LJs  - cage_in_solvent_LJ

list_E_bound_C.append(E_bound_C)
list_E_bound_LJ.append(E_bound_LJ)

print("list_E_bound_C= ", list_E_bound_C)
print("list_E_bound_LJ= ", list_E_bound_LJ)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G3/free/discharge/output/lambda-0.00')

print("---------------- Cage G3 ----------------")


hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul

hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs
sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()


E_free_C = total_host_guest_Coul + sol_sol_Coul - sol_sol_Coul
E_free_LJ = total_host_guest_LJs + sol_sol_LJs - sol_sol_LJs

list_E_free_C.append(E_free_C)
list_E_free_LJ.append(E_free_LJ)

print("list_E_free_C= ", list_E_free_C)
print("list_E_free_LJ= ", list_E_free_LJ)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G3/bound/discharge/output/lambda-0.00')

hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul
hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs

sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()

E_bound_C = total_host_guest_Coul  - cage_in_solvent_elec
E_bound_LJ = total_host_guest_LJs  - cage_in_solvent_LJ

list_E_bound_C.append(E_bound_C)
list_E_bound_LJ.append(E_bound_LJ)

print("list_E_bound_C= ", list_E_bound_C)
print("list_E_bound_LJ= ", list_E_bound_LJ)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G4/free/discharge/output/lambda-0.00')

print("---------------- Cage G4 ----------------")


hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul

hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs
sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()


E_free_C = total_host_guest_Coul + sol_sol_Coul - sol_sol_Coul
E_free_LJ = total_host_guest_LJs + sol_sol_LJs - sol_sol_LJs

list_E_free_C.append(E_free_C)
list_E_free_LJ.append(E_free_LJ)

print("list_E_free_C= ", list_E_free_C)
print("list_E_free_LJ= ", list_E_free_LJ)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G4/bound/discharge/output/lambda-0.00')

hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul
hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs

sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()

E_bound_C = total_host_guest_Coul  - cage_in_solvent_elec
E_bound_LJ = total_host_guest_LJs  - cage_in_solvent_LJ

list_E_bound_C.append(E_bound_C)
list_E_bound_LJ.append(E_bound_LJ)

print("list_E_bound_C= ", list_E_bound_C)
print("list_E_bound_LJ= ", list_E_bound_LJ)


print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G5/macromodel/mod_NH2/test1/free/discharge/output/lambda-0.0')

print("---------------- Cage G5 ----------------")


hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul

hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs
sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()


E_free_C = total_host_guest_Coul + sol_sol_Coul - sol_sol_Coul
E_free_LJ = total_host_guest_LJs + sol_sol_LJs - sol_sol_LJs

list_E_free_C.append(E_free_C)
list_E_free_LJ.append(E_free_LJ)

print("list_E_free_C= ", list_E_free_C)
print("list_E_free_LJ= ", list_E_free_LJ)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G5/macromodel/mod_NH2/test1/bound/discharge/output/lambda-0.0')


hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul
hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs

sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()

E_bound_C = total_host_guest_Coul  - cage_in_solvent_elec
E_bound_LJ = total_host_guest_LJs  - cage_in_solvent_LJ

list_E_bound_C.append(E_bound_C)
list_E_bound_LJ.append(E_bound_LJ)

print("list_E_bound_C= ", list_E_bound_C)
print("list_E_bound_LJ= ", list_E_bound_LJ)


print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G6/free/discharge/output/lambda-0.00')

print("---------------- Cage G6 ----------------")


hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul

hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs
sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()


E_free_C = total_host_guest_Coul + sol_sol_Coul - sol_sol_Coul
E_free_LJ = total_host_guest_LJs + sol_sol_LJs - sol_sol_LJs

list_E_free_C.append(E_free_C)
list_E_free_LJ.append(E_free_LJ)

print("list_E_free_C= ", list_E_free_C)
print("list_E_free_LJ= ", list_E_free_LJ)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G6/bound/discharge/output/lambda-0.00')

hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul
hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs

sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()

E_bound_C = total_host_guest_Coul  - cage_in_solvent_elec
E_bound_LJ = total_host_guest_LJs  - cage_in_solvent_LJ

list_E_bound_C.append(E_bound_C)
list_E_bound_LJ.append(E_bound_LJ)

print("list_E_bound_C= ", list_E_bound_C)
print("list_E_bound_LJ= ", list_E_bound_LJ)



print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G7/free/discharge/output/lambda-0.00')

print("---------------- Cage G7 ----------------")


hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul

hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs
sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()


E_free_C = total_host_guest_Coul + sol_sol_Coul - sol_sol_Coul
E_free_LJ = total_host_guest_LJs + sol_sol_LJs - sol_sol_LJs

list_E_free_C.append(E_free_C)
list_E_free_LJ.append(E_free_LJ)

print("list_E_free_C= ", list_E_free_C)
print("list_E_free_LJ= ", list_E_free_LJ)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G7/bound/discharge/output/lambda-0.00')

hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul
hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs

sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()

E_bound_C = total_host_guest_Coul  - cage_in_solvent_elec
E_bound_LJ = total_host_guest_LJs  - cage_in_solvent_LJ

list_E_bound_C.append(E_bound_C)
list_E_bound_LJ.append(E_bound_LJ)

print("list_E_bound_C= ", list_E_bound_C)
print("list_E_bound_LJ= ", list_E_bound_LJ)


print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G8/free/discharge/output/lambda-0.00')

print("---------------- Cage G8 ----------------")


hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul

hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs
sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()


E_free_C = total_host_guest_Coul + sol_sol_Coul - sol_sol_Coul
E_free_LJ = total_host_guest_LJs + sol_sol_LJs - sol_sol_LJs

list_E_free_C.append(E_free_C)
list_E_free_LJ.append(E_free_LJ)

print("list_E_free_C= ", list_E_free_C)
print("list_E_free_LJ= ", list_E_free_LJ)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

os.chdir( '/media/sofia/LACIE SHARE/3rd_year/final_cages_files/plato_FEP/cage_G8/bound/discharge/output/lambda-0.00')

hard_sol_Coul =hard_sol_C()
hard_intraCoul = hard_intraC()
hard_intraff14Coul = hard_intraff14C()
total_host_guest_Coul = hard_sol_Coul + hard_intraCoul + hard_intraff14Coul
hard_sol_LJs =hard_sol_LJ()
hard_intraLJs = hard_intraLJ()
hard_intraff14LJs = hard_intraff14LJ()
total_host_guest_LJs = hard_sol_LJs + hard_intraLJs + hard_intraff14LJs

sol_sol_LJs =sol_sol_LJ()
sol_sol_Coul = sol_sol_C()

E_bound_C = total_host_guest_Coul  - cage_in_solvent_elec
E_bound_LJ = total_host_guest_LJs  - cage_in_solvent_LJ

list_E_bound_C.append(E_bound_C)
list_E_bound_LJ.append(E_bound_LJ)

print("list_E_bound_C= ", list_E_bound_C)
print("list_E_bound_LJ= ", list_E_bound_LJ)

diff_C = []
diff_LJ =[]
for i in range(0,len(list_E_bound_C)):
    diff_C.append(list_E_bound_C[i] - list_E_free_C[i])
    diff_LJ.append(list_E_bound_LJ[i] - list_E_free_LJ[i])


exper = [-7.51, -5.23, -10.41 ,-12.06 ,-10.41, -5.40 ,-3.62, -3.46]
ddg_mod_char2 = [-11.95 ,-13.78 ,-5.90 ,-2.64 ,-1.41 ,-16.02 ,-9.79 ,-6.48]
target = pd.DataFrame(ddg_mod_char2, columns = ["DDG"])
experiment = pd.DataFrame(exper, columns = ["DDG"])
# df = pd.DataFrame({"elec": diff_C, "LJ": diff_LJ })
df = pd.DataFrame({"elec_bound": list_E_bound_C, "LJ_bound": list_E_bound_LJ, "elec_free": list_E_free_C, "LJ_free": list_E_free_LJ})

X = df[["elec_bound", "LJ_bound","elec_free", "LJ_free" ]]
# X = df[["elec", "LJ"]]
y = experiment["DDG"]
X = sm.add_constant(X)
# model = sm.OLS(y, X).fit()
# predictions = model.predict(X) # make the predictions by the model
# Print out the statistics

model = sm.OLS(y, X).fit()
predictions = model.predict(X)

print(model.summary())
print(model.params)
# print(model.bse) #std error

constant = model.params[0]
elec_b_coef = model.params[1]
LJ_b_coef = model.params[2]
elec_f_coef = model.params[3]
LJ_f_coef = model.params[4]

predicted_bound_DDG = []
for i in range (0, len(list_E_bound_LJ)):
    # pred = elec_b_coef*diff_C[i] + LJ_b_coef* diff_LJ[i] +constant
    pred = elec_b_coef*list_E_bound_C[i] + LJ_b_coef* list_E_bound_LJ[i]  +elec_f_coef*list_E_free_C[i] + LJ_f_coef* list_E_free_LJ[i] +constant
    predicted_bound_DDG.append(pred)

df_res_bound = pd.DataFrame({"Experimenal": ddg_mod_char2, "Predicted": predicted_bound_DDG})
print(df_res_bound )
