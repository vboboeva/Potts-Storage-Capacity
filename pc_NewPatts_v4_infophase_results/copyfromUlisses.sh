#!/bin/bash



for i in 70; #1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900;
  do
    dir='zeta0';
    a=`echo $i \* 0.001 |bc`;
    echo $dir$a;
	scp vboboeva@frontend1.hpc.sissa.it:/scratch/vboboeva/pc_corrv4/$dir$a/storage_capacity_dataset_0 data/storage_capacity_$dir$a
  done
  
#for i in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 500;
  #do
    #dir='zeta';
    #a=`echo $i \* 1.000 |bc`;
    #echo $dir$a;
	#scp vboboeva@frontend1.hpc.sissa.it:/scratch/vboboeva/pc_corrv4/$dir$a/storage_capacity_dataset_0 data/storage_capacity_$dir$a
  #done  
