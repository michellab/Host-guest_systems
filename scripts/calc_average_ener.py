

def hard_sol_C(): 
    hard_solv_C = []
    count = 0
    iname=("energies.txt")
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
    iname=("energies.txt")
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
    iname=("energies.txt")
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
    iname=("energies.txt")
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
    iname=("energies.txt")
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
    iname=("energies.txt")
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



hard_sol_C =hard_sol_C()
hard_intraC = hard_intraC()
hard_intraff14C = hard_intraff14C()
total_host_guest_C = hard_sol_C + hard_intraC + hard_intraff14C
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("solute- solvent electrostatic interactions:", hard_sol_C)3
print("intramolecular electrostatic interactions:", hard_intraC)
print("intramolecular 1-4 electrostatic interactions:", hard_intraff14C)
print("Total electrostatic interactions:", total_host_guest_C)
print("---------------------------------")
print("---------------------------------")

hard_sol_LJ =hard_sol_LJ()
hard_intraLJ = hard_intraLJ()
hard_intraff14LJ = hard_intraff14LJ()
total_host_guest_LJ = hard_sol_LJ + hard_intraLJ + hard_intraff14LJ
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("solute- solvent LJ interactions:", hard_sol_LJ)
print("intramolecular LJ interactions:", hard_intraLJ)
print("intramolecular 1-4 LJ interactions:", hard_intraff14LJ)
print("Total LJ interactions:", total_host_guest_LJ)


def sol_sol_LJ(): 
    solv_solv_LJ = []
    count = 0
    iname=("energies.txt")
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
    iname=("energies.txt")
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



def sol_intra_LJ(): 
    solv_intra_LJ = []
    count = 0
    iname=("energies.txt")
    ener_list = []
    with open(iname) as infile:
        lines = infile.readlines()
        for line in lines:
            for i in range (0, len(line.split()) ):
                if line.split()[i] == "E_{solvent_intraclj}^{LJ}":
                    solv_intra_LJ.append(line.split()[i+2].strip(","))
                    entry = float(line.split()[i+2].strip(","))
                    count = count + entry
    # print ("solv_intra_LJ = ",solv_intra_LJ)
    avg = count/len(solv_intra_LJ)
    print("Average E_{solvent_intraclj}^{LJ} = ", avg)
    return(avg)

       
def sol_intra_C(): 
    solv_intra_C = []
    count = 0
    iname=("energies.txt")
    ener_list = []
    with open(iname) as infile:
        lines = infile.readlines()
        for line in lines: 
            for i in range (0, len(line.split()) ):
                if line.split()[i] == "E_{solvent_intraclj}^{coulomb}":
                    solv_intra_C.append(line.split()[i+2].strip(","))
                    entry = float(line.split()[i+2].strip(","))
                    count = count + entry
    # print ("solv_intra_C = ",solv_intra_C)
    avg = count/len(solv_intra_C)
    print("Average E_{solvent_intraclj}^{coulomb} = ", avg)
    return(avg)


sol_sol_LJ =sol_sol_LJ()
sol_sol_C = sol_sol_C()
sol_intra_C = sol_intra_C()
sol_intra_LJ = sol_intra_LJ()
print()
print("~~~~~~~~~~~ SOLVENT JUST FOR FREE STATE ~~~~~~~~~~~~~~~~~~~~")
print("Average E_{solvent:solvent}^{coulomb}:", sol_sol_C)
print("Average E_{solvent:solvent}^{LJ}:", sol_sol_LJ)
print("Average E_{solvent_intraclj}^{coulomb}:", sol_intra_C)
print("Average E_{solvent_intraclj}^{LJ}:", sol_intra_LJ)

