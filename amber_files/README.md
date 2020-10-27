System combination with BSS:
```
run ./combine.py --system1 ../G4/g4q.prm7 ../G4/g4q.rst7 --system2 ../c2v_ligpargen/lig_host1_from_solv.prm7 ../c2v_ligpargen/lig_host1_from_solv.rst7 --output guest_1h

run ./combine.py --system1 guest_1h.prm7 guest_1h.rst7 --system2 ../c2v_ligpargen/lig_host2_from_solv.prm7 ../c2v_ligpargen/lig_host2_from_solv.rst7 --output guest_2h

run ./combine.py --system1 guest_1h.prm7 guest_1h.rst7 --system2 ../c2v_ligpargen/lig_host2_from_solv.prm7 ../c2v_ligpargen/lig_host2_from_solv.rst7 --output guest_2                                                                                             

run ./combine.py --system1 guest_2.prm7 guest_2.rst7 --system2 ../c2v_ligpargen/lig_host3_from_solv.prm7 ../c2v_ligpargen/lig_host3_from_solv.rst7 --output guest_3            

run ./combine.py --system1 guest_3.prm7 guest_3.rst7 --system2 ../c2v_ligpargen/lig_host4_from_solv.prm7 ../c2v_ligpargen/lig_host4_from_solv.rst7 --output guest_4            

run ./combine.py --system1 guest_4.prm7 guest_4.rst7 --system2 ../c2v_ligpargen/mtl_from_solv1.prm7 ../c2v_ligpargen/mtl_from_solv1.rst7 --output guest_host_mtl    

run ./combine.py --system1 guest_host_mtl.prm7 guest_host_mtl.rst7 --system2 ../c2v_ligpargen/mtl_from_solv2.prm7 ../c2v_ligpargen/mtl_from_solv2.rst7 --output guest_host_mtl2

run ./combine.py --system1 guest_host_mtl2.prm7 guest_host_mtl2.rst7 --system2 ../barf/barf1_from_solvated.prm7 ../barf/barf1_from_solvated.rst7 --output guest_host_mtl_b1

run ./combine.py --system1 guest_host_mtl_b1.prm7 guest_host_mtl_b1.rst7 --system2 ../barf/barf2_from_solvated.prm7 ../barf/barf2_from_solvated.rst7 --output guest_host_mtl_b2

run ./combine.py --system1 guest_host_mtl_b2.prm7 guest_host_mtl_b2.rst7 --system2 ../barf/barf3_from_solvated.prm7 ../barf/barf3_from_solvated.rst7 --output guest_host_mtl_b3

run ./combine.py --system1 guest_host_mtl_b3.prm7 guest_host_mtl_b3.rst7 --system2 ../barf/barf4_from_solvated.prm7 ../barf/barf4_from_solvated.rst7 --output guest_host_mtl_b4

run ./combine.py --system1 guest_host_mtl_b4.prm7 guest_host_mtl_b4.rst7 --system2 DCM_from_solv.prm7 DCM_from_solv.rst7 --output cage_G4_solv
 ```
 
