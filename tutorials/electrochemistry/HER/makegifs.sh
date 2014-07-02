#!/bin/bash

for i in rate coverage production_rate
	do convert -delay 50 -loop 0 voltage0.0/${i}0.0.png voltage-*/${i}*.png ${i}.gif
done