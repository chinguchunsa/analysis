V33 :0x4 utils
13 mod_utils.f90 S624 0
12/09/2019  13:23:37
enduse
D 56 21 9 2 12 16 0 0 1 0 0
 0 14 3 3 15 15
 0 14 15 3 15 15
D 59 21 12 1 3 19 0 0 1 0 0
 0 18 3 3 19 19
D 62 21 12 1 3 15 0 0 1 0 0
 0 14 3 3 15 15
D 65 21 7 2 20 25 0 0 1 0 0
 0 22 3 3 23 23
 0 24 23 3 24 24
D 68 21 7 2 20 26 0 0 1 0 0
 0 22 3 3 23 23
 0 3 23 3 3 3
D 71 21 12 1 3 23 0 0 1 0 0
 0 22 3 3 23 23
D 74 21 9 2 27 31 0 0 1 0 0
 0 29 3 3 30 30
 0 29 30 3 30 30
D 77 21 9 2 27 31 0 0 1 0 0
 0 29 3 3 30 30
 0 29 30 3 30 30
S 624 24 0 0 0 8 1 0 5015 10005 0 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 utils
S 625 23 5 0 0 0 626 624 5021 0 0 A 0 0 0 0 B 0 46 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_indices
S 626 14 5 0 0 0 1 625 5021 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 10 0 624 0 0 0 0 get_indices
F 626 0
S 627 23 5 0 0 0 628 624 5033 0 0 A 0 0 0 0 B 0 106 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_indices_ip
S 628 14 5 0 0 0 1 627 5033 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 49 0 624 0 0 0 0 get_indices_ip
F 628 0
S 629 23 5 0 0 0 636 624 5048 0 0 A 0 0 0 0 B 0 130 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_psi_det
S 630 1 3 3 0 7 1 629 5060 4 3000 A 0 0 0 0 B 0 130 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 i
S 631 6 3 1 0 7 1 629 5062 800004 3000 A 0 0 0 0 B 0 130 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nstuse
S 632 6 3 1 0 7 1 629 5069 800004 3000 A 0 0 0 0 B 0 130 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nstates
S 633 7 3 1 0 56 1 629 5077 800204 3000 A 0 0 0 0 B 0 130 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ci_vec
S 634 7 3 1 0 59 1 629 5084 800204 3000 A 0 0 0 0 B 0 130 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 psi_itime
S 635 7 3 3 0 62 1 629 5094 800204 3000 A 0 0 0 0 B 0 130 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 psi_det
S 636 14 5 0 0 0 1 629 5048 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4 6 0 0 0 0 0 0 0 0 0 0 0 0 109 0 624 0 0 0 0 get_psi_det
F 636 6 630 631 632 633 634 635
S 637 6 1 0 0 6 1 629 5102 40800006 3000 A 0 0 0 0 B 0 130 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_14
S 638 6 1 0 0 6 1 629 5109 40800006 3000 A 0 0 0 0 B 0 130 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_16
S 639 6 1 0 0 6 1 629 5116 40800006 3000 A 0 0 0 0 B 0 130 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_18
S 640 6 1 0 0 6 1 629 5123 40800006 3000 A 0 0 0 0 B 0 130 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_21
S 641 23 5 0 0 0 651 624 5130 0 0 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_1density_ip
S 642 7 3 1 0 71 1 641 5094 800204 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 psi_det
S 643 6 3 0 0 7 1 641 5069 800004 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nstates
S 644 7 3 1 0 65 1 641 5146 800204 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 hole
S 645 7 3 1 0 68 1 641 5151 800204 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 part
S 646 7 3 3 0 74 1 641 5156 800204 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rho_a
S 647 7 3 3 0 77 1 641 5162 800204 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rho_b
S 648 6 3 0 0 7 1 641 5168 800004 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nrorb
S 649 1 3 0 0 7 1 641 5174 4 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 noa
S 650 1 3 0 0 7 1 641 5178 4 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nob
S 651 14 5 0 0 0 1 641 5130 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 11 9 0 0 0 0 0 0 0 0 0 0 0 0 133 0 624 0 0 0 0 get_1density_ip
F 651 9 642 643 644 645 646 647 648 649 650
S 652 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 653 6 1 0 0 6 1 641 5182 40800006 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_22
S 654 6 1 0 0 6 1 641 5189 40800006 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_25
S 655 6 1 0 0 6 1 641 5196 40800006 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_27
S 656 6 1 0 0 6 1 641 5203 40800006 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_24
S 657 6 1 0 0 6 1 641 5210 40800006 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_31
S 658 6 1 0 0 6 1 641 5217 40800006 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_33
S 659 6 1 0 0 6 1 641 5224 40800006 3000 A 0 0 0 0 B 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_35
A 12 1 0 0 0 6 639 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 13 1 0 0 0 7 632 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 14 7 0 0 0 6 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15 1 0 0 0 6 637 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 16 1 0 0 0 6 638 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 17 1 0 0 0 7 631 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 18 7 0 0 0 6 17 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 19 1 0 0 0 6 640 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 20 1 0 0 0 6 655 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 21 1 0 0 0 7 643 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 22 7 0 0 0 6 21 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 23 1 0 0 0 6 653 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 24 2 0 0 0 6 652 0 0 0 24 0 0 0 0 0 0 0 0 0 0 0
A 25 1 0 0 0 6 654 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 26 1 0 0 0 6 656 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 27 1 0 0 0 6 659 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 28 1 0 0 0 7 648 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 29 7 0 0 0 6 28 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 30 1 0 0 0 6 657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 31 1 0 0 0 6 658 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
Z
