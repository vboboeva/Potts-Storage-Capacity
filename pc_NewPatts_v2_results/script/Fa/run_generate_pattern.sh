 #!/bin/bash
 
for i in  `seq 1 9`; 
	do
		a=`echo $i \* 0.1 |bc`;
		mpic++ generate_pattern.cpp -o generate_pattern;
		mpirun -n 1 generate_pattern $a;
	done
	

