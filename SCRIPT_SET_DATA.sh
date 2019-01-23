#!/bin/sh

NUMBER_OF_THREADS_TO_BE_USED="$1"	# number of processors

INPUT_1="$2"								# Do you want to random select a new instance?
INPUT_2="$3"								# Compute the random seed and rewrite the files
												# used in the search of smooth divisors
																								
echo ""
echo "$NUMBER_OF_THREADS_TO_BE_USED" > ./inc/NUMBER_OF_THREADS
echo "The file ./inc/NUMBER_OF_THREADS has been updated!";

if [ "$INPUT_1" = "NEW_INSTANCE" ] || [ "$INPUT_2" = "NEW_INSTANCE" ]; then
	(magma main/GET_INSTANCE.m) > ./outputs/OUTPUT_OBTAINED_BY_THE_SCRIPT_GET_INSTANCE
	echo "NEW RANDOM INSTANCE GENERATED!"
else
	echo "OLD RANDOM INSTANCE READY TO BE USED!"
fi

if [ "$INPUT_1" = "NEW_RANDOM_WALK_SEED" ] || [ "$INPUT_2" = "NEW_RANDOM_WALK_SEED" ]; then
	(magma main/GET_SEED.m) > ./outputs/OUTPUT_OBTAINED_BY_THE_SCRIPT_GET_SEED
	echo "The files ./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_* have been updated..."
	echo "NEW RANDOM WALK GENERATED!"
else
	echo "OLD RANDOM WALK READY TO BE USED!"
fi

(magma main/GET_MATRIX.m) > ./outputs/OUTPUT_OBTAINED_BY_THE_SCRIPT_GET_MATRIX
echo "The files ./outputs/GET_MATRIX_POSITION_USING_THREAD_* have been updated..."
echo "";

echo "ALL THE DATA REQUIRED IN THE PARALLEL VERSION OF THE GAUDRY ALGORITHM (INDEX CALCULUS ALGORITHM) IS READY."
echo "You should run the following commands:"
echo "   sh SCRIPT_GET_SMOOTH_DIVISORS.sh"
echo "   sh SCRIPT_GET_MATRIX_POSITION.sh"
echo "   sh SCRIPT_JOIN_THE_MATRIX.sh"
echo "   sh SCRIPT_MATRIX_TRANSPOSE.sh"
echo "   sh SCRIPT_VECTOR_SOLUTION.sh"
echo ""
