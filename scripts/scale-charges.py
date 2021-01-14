
import sys

def scale_charges(scale_factor, num_atoms):
    scale_factor = float(scale_factor)
    num_atoms = int(num_atoms)
    print(type(scale_factor))
    print(type(num_atoms))

    iname=("SYSTEM_bck.top")
    ch_list = []
    with open(iname) as infile:
        while True:
            line=infile.readline()
            if line=="":
                break
            elif "%FLAG CHARGE" in line:
                line=infile.readline()
                while "%FLAG ATOMIC_NUMBER" not in line:
                    line=infile.readline()
                    ch_list.append(line.split())

    ch_list= ch_list[:-1]

    charges_list=[]
    for i in ch_list:
        for j in i: 
            charges_list.append(float(j))

  
    new_charges=[]
    count = 0
    while count <= num_atoms-1:
        new_c = charges_list[count]/scale_factor
        new_charges.append(new_c)
        count = count+1

    for char in charges_list[num_atoms :]:
        new_charges.append(char)

    index =[]
    for i in range(0, len(new_charges)):
        index.append(i)

    with open(iname) as infile:
        
        with open("SYSTEM_mod.top", "w") as o:
        
            while True:
                line=infile.readline()
                if line=="":
                    break
                elif "%FLAG CHARGE" not in line:
                    o.write(line)
                else:
                    o.write(line)
                    line=infile.readline()
                    o.write(line)
                    #we need to repeat the following code 4 times to get all the charges with the correct format
                    count = 0
                    while count < len(new_charges):
                        if count < len(new_charges) and new_charges[count] < 0:
                            o.write(' '+str(format(new_charges[count], ".8E")))
                            count = count+1
                            print(count)
                        elif count < len(new_charges) and new_charges[count] >=0:
                            o.write('  '+str(format(new_charges[count], ".8E")))
                            count = count+1
                            print(count)
                        elif  count > len(new_charges):
                            break

                        if count < len(new_charges) and new_charges[count] < 0:
                            o.write(' '+str(format(new_charges[count], ".8E")))
                            count = count+1
                            print(count)
                        elif count < len(new_charges) and new_charges[count] >=0:
                            o.write('  '+str(format(new_charges[count], ".8E")))
                            count = count+1
                            print(count)
                        elif  count > len(new_charges):
                            break
                        
                        if count < len(new_charges) and new_charges[count] < 0:
                            o.write(' '+str(format(new_charges[count], ".8E")))
                            count = count+1
                            print(count)
                        elif count < len(new_charges) and new_charges[count] >=0:
                            o.write('  '+str(format(new_charges[count], ".8E")))
                            count = count+1
                            print(count)
                        elif  count > len(new_charges):
                            break
                        
                        if count < len(new_charges) and new_charges[count] < 0:
                            o.write(' '+str(format(new_charges[count], ".8E")))
                            count = count+1
                            print(count)
                        elif count < len(new_charges) and new_charges[count] >=0:
                            o.write('  '+str(format(new_charges[count], ".8E")))
                            count = count+1
                            print(count)
                        elif  count > len(new_charges):
                            break

                        if count < len(new_charges)  and new_charges[count] < 0:
                            o.write(' '+str(format(new_charges[count], ".8E"))+'\n')
                            count = count+1
                            print(count)
                        elif count < len(new_charges)  and new_charges[count] >=0:
                            o.write('  '+str(format(new_charges[count], ".8E"))+'\n')
                            count = count+1
                            print(count)
                        elif  count > len(new_charges):
                            break
                   
                    o.write('\n')
                    while "%FLAG ATOMIC_NUMBER" not in line:
                        line=infile.readline()
                    o.write(line)


if __name__ == '__main__':
    scale_factor = sys.argv[1]
    num_atoms = sys.argv[2]
  
    scale_charges(scale_factor, num_atoms)

