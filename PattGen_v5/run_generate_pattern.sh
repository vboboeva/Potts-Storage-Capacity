 #!/bin/bash
g++ generate_pattern.cpp -o generate_pattern -lm
./generate_pattern

#mpic++ compute_corrs_patts.cpp -o compute_corrs_patts
#mpirun -n 1 compute_corrs_patts 

#mpic++ compute_corrs_units.cpp -o compute_corrs_units
#mpirun -n 1 compute_corrs_units 

	

