#!/bin/bash

cd zeta0.000/ && bash compute_pc.sh && cd ../

for z in  1 2 5 8 9 10 11 12 15 20 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950;
  do
    dir='zeta0';
    a=`echo $z \* 0.001 |bc`;
    cd $dir$a/;
    bash compute_pc.sh;
    cd ../
  done
  
  


