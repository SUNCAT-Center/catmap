#!/bin/bash

mkdir $1
# for i in -1.3 -1.1 -0.9 -0.7 -0.5 -0.3 -0.1 0.1 0.3 0.5 0.7 0.9 1.1 1.3
for i in 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.76 0.77 0.78 0.79 0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.0
 do 
 if [ -d "$1/voltage$i" ]; then
 	echo $i "is done"
 	continue
 fi
 mkdir $1/voltage$i
 cd $1/voltage$i
 cp ../../mkm_job.py .
 cp ../../ORR_template.mkm .
 cp ../../make_input.py .
 python make_input.py
 python mkm_job.py $i
 cd ../..
done
