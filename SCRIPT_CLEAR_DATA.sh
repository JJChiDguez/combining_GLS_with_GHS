#!/bin/sh

echo ""
INPUT_1="$1"			# Do you want to delete the data about the instance?
INPUT_2="$2"			# Do you want to delete the data about the Random Walk?

# DELETE ALL THE  PREVIOUS DATA ABOUT THE MATRIX AND THE VECTOR IN THE KERNEL
if [ "$INPUT_1" = "INSTANCE" ] || [ "$INPUT_2" = "INSTANCE" ]; then
	echo "The folder ./inc/tmp is empty now!"
	rm -f ./inc/tmp/*
fi;

echo "The folder ./outputs/tmp is empty now!"
rm -f ./outputs/tmp/*

# DELETE THE RANDOM SEED AND THE FILES USED TO COMPUTE PREVIOUS DATA ABOUT 
# THE MATRIX AND THE VECTOR IN THE KERNEL
rm -f ./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_*
rm -f ./outputs/GET_MATRIX_POSITION_USING_THREAD_*

if [ "$INPUT_1" = "RANDOM_WALK_SEED" ] || [ "$INPUT_2" = "RANDOM_WALK_SEED" ]; then
	echo "The files ./outputs/RANDOM_WALK and ./outputs/T0 have been deleted!"
	rm -f ./outputs/RANDOM_WALK
	rm -f ./outputs/T0
fi;

# DELETE ALL PREVIOUS LOGS
echo "The folder ./outputs/logs is empty now!"
rm -f ./outputs/logs/*

# DELETE THE SCRIPTS USED IN PREVIOUS RUNS
rm -f SCRIPT_GET_SMOOTH_DIVISORS.sh
rm -f SCRIPT_GET_MATRIX_POSITION.sh

# Logs for the parallel run
rm -f ./outputs/OUTPUT_OBTAINED_BY_THE_SCRIPT_GET_*
echo ""
