#!/bin/bash



for i in 1100 1500 2000 3000 4000 5000 10000; #1 2 5 8 9 10 11 12 15 20 50 100 150 200 500 800; #
  do
    z=`echo $i \* 0.001000 |bc`;
	echo $z;
	mkdir zeta$z;
	scp vboboeva@frontend1.hpc.sissa.it:/scratch/vboboeva/Reviews/zeta$z/storage_capacity_dataset_0  zeta$z/storage_capacity_dataset_0
  done
