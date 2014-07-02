#!/bin/bash
#
mkdir $1
# for i in -1.3 -1.1 -0.9 -0.7 -0.5 -0.3 -0.1 0.1 0.3 0.5 0.7 0.9 1.1 1.3
for i in 0.0 -0.05 -0.1 -0.15 -0.2 -0.25 -0.3 -0.35 -0.4 -0.45 -0.5 -0.55 -0.6 -0.65 -0.7 -0.75 -0.8 -0.85 -0.9 -0.95 -1.0
 do mkdir $1/voltage$i
 cd $1/voltage$i
 cp ../../do_kinetics.py .
 cp ../../CO2_reduction_template.mkm .
 cp ../../make_input.py .
 python make_input.py $i
 python2.7 do_kinetics.py $i
 cd ../..
done
# cat CO2_reduction_template.mkm | sed "s/voltage \= -0.5/voltage \= $1/" > CO2_reduction.mkm
# python2.7 red_job.py $1
