-------------------------------------------------------------------------------

Extending the GLS endomorphism to speedup the GHS Weil descent using Magma


This research project was realized by:

1. Chi-Domínguez Jesús-Javier,
2. Rodríguez-Henríquez Francisco, and 
3. Benjamin Smith

-------------------------------------------------------------------------------

Code created by Jesús Javier Chi Domínguez, <jjchi@computacion.cs.cinvestav.mx>,
											and <chidoys@gmail.com>

-------------------------------------------------------------------------------


CONTENTS OF THIS FILE

---------------------

	* Introduction
	* Requirements
	* Remarks
	* Configuration
 
---------------------


-------------------------------------------------------------------------------

INTRODUCTION

This a Magma code implementation of the Enge-Gaudry algorithm was implemented 
(an Index-Calculus based algorithm), and an efficient endomorphism on the 
jacobian of the image of the GHS Weil descent technique. In addition, our 
Index-Calculus based algorithm uses that endomorphism in order to speedup the 
Discrete Logarithm computation. All the implementation is over finite field 
extensions of GF(2). 



-------------------------------------------------------------------------------

REQUIREMENTS

This program requires to have installed the algebra system Magma.



-------------------------------------------------------------------------------

REMARKS

Before running our Magma code implementation, one should create the folders 
outputs/tmp and outputs/logs. This step is done only one time by running the 
following two commands:

	mkdir inc/tmp
	mkdir outputs/tmp
	mkdir outputs/logs



-------------------------------------------------------------------------------

CONFIGURATION

[FINDING DISCRETE LOGARITHMS: GLS combined with GHS]

First, one needs to create the files that contain the following:

	--------------------------------------------------------------------------
	1.	the finite field to be used which must be described as in the files 
		in folder ./inc/GF.
			THE FILE ./inc/FINITE_FIELD.m MUST BE UPDATED,
	2.	the elliptic curve to be used which must be described as in the files 
		in folder ./inc/EC. In addition, 
			THE FILE ./inc/ELLIPTIC_CURVE.m must MUST BE UPDATED,
	3.	the points P and P' such that P' = [\lambda] P,
	4.	the genus-g hyperelliptic curve that the GHS Weil descent technique 
		constructs,
	5.	the divisors D = GHS(P) and D' = GHS(P')

	TO CREATE OR UPDATE THOSE PREVIUS DATA, TO RUN: magma main/GET_INSTANCE.m
	--------------------------------------------------------------------------
	6.	the number of threads to be used in the search of smooth divisors must 
		be set in the file ./inc/NUMBER_OF_THREADS,
	7.	the number of threads to be used in the linear algebra step, smoothness 
		bound, and number of random divisors in the sequence (i.e., the list 
		to be used for the random walk).
			THE FILE ./inc/SMOOTHNESS_PARAMETERS.m MUST BE UPDATED, and
	8.	The random walk seed required in the Enge-Gaudry algorithm.
	
	TO CREATE OR UPDATE THOSE PREVIUS DATA, TO RUN:	magma main/GET_SEED.m

	However, you can update all the previous data by running the following
		sh SCRIPT_SET_DATA.sh NUMBER_OF_THREADS INPUT_1 INPUT_2,
	where INPUT_1 and INPUT_2 can take the values NEW_RANDOM_WALK_SEED or 
	NEW_INSTANCE. And then, by modifying the file 
		./inc/SMOOTHNESS_PARAMETERS.m
	
Once that one has the framework correctly set. The random instances can be solved 
by running the following
	
	sh SCRIPT_GET_SMOOTH_DIVISORS.sh
	sh SCRIPT_GET_MATRIX_POSITION.sh
	sh SCRIPT_JOIN_THE_MATRIX.sh
	sh SCRIPT_MATRIX_TRANSPOSE.sh
	sh SCRIPT_VECTOR_SOLUTION.sh
		
In addition, to delete the old data you should run the following

	sh SCRIPT_CLEAR_DATA.sh INPUT_3 INPUT_4
		
where INPUT_3 and INPUT_4 can take the values RANDOM_WALK_SEED or INSTANCE


REMARKS: The following info could save your day:
	
	If one of your process created by script SCRIPT_GET_SMOOTH_DIVISORS.sh die
	and you want to use the DATA saved, then you can run the following script:
		sh SCRIPT_CHECK_SIZES.sh THREAD_ID
	This script computes the number of s-smooth divisors found before the i-th
	thread died


[TESTING THE GENERALIZED GLS ENDOMORPHISM IN THE JACOBIAN]


One needs to create the files that contain the following:	

	1.	the finite field to be used which must be described as in the files 
		in folder ./inc/GF.
			THE FILE ./inc/FINITE_FIELD.m MUST BE UPDATED,
	2.	the elliptic curve to be used which must be described as in the files 
		in folder ./inc/EC. In addition, 
			THE FILE ./inc/ELLIPTIC_CURVE.m must MUST BE UPDATED,
	3.	the points P and P' such that P' = [\lambda] P,
	4.	the genus-g hyperelliptic curve that the GHS Weil descent technique 
		constructs,
	5.	the divisors D = GHS(P) and D' = GHS(P')

	TO CREATE OR UPDATE THOSE PREVIUS DATA, TO RUN: magma main/GET_INSTANCE.m
	
Once that one has the framework correctly set, to run:
	magma lib/ENDOMORPHISM_IN_THE_JACOBIAN.m
