/* ----------------------------------------------------------------------------- */
load "./inc/FINITE_FIELD.m";
/* ----------------------------------------------------------------------------- */
load "./inc/tmp/ELLIPTIC_CURVE_parameters";
E_qn := EllipticCurve([F_qn| 1, a_qn, 0, 0, b_qn]);
E_qn_size := c*r;
load "./inc/tmp/ELLIPTIC_CURVE_instance";
load "./inc/SMOOTHNESS_PARAMETERS.m";

/* ----------------------------------------------------------------------------- */
load "./inc/tmp/GENUS_g_HCURVE_parameters";
H_q := HyperellipticCurve(f_q, h_q);
J_q := Jacobian(H_q);
load "./inc/tmp/GENUS_g_HCURVE_instance";

/* ----------------------------------------------------------------------------- */
load "./lib/LIBRARY.m";
I_time := Cputime();

System("wc -l < ./outputs/tmp/INDXS > ./outputs/tmp/N_size");
N := StringToInteger(Gets(Open("./outputs/tmp/N_size", "r")));
M := SparseMatrix(N, N);

System("wc -l < ./outputs/tmp/tr_INDXS > ./outputs/tmp/L_size");
L := StringToInteger(Gets(Open("./outputs/tmp/L_size", "r")));

vlstr_lst := GetingSmoothnessCoefFromStringFile("./outputs/tmp/tr_COEFS", L);
idxtr_lst := GetingSmoothnessCoefFromStringFile("./outputs/tmp/tr_INDXS", L);


for var_j in [1 .. L] do
	for var_k in [1 .. #vlstr_lst[var_j]] do
		M[var_j, idxtr_lst[var_j, var_k]] +:= vlstr_lst[var_j, var_k];
	end for;
end for;

ALP_lst := GetingSmoothnessCoefFromStringFile("./outputs/tmp/ALP", N);
BTA_lst := GetingSmoothnessCoefFromStringFile("./outputs/tmp/BTA", N);


F_time := Cputime();
printf "> %o s,\nWe proceed to find a vector of the kernel of an \n%o\t", F_time - I_time, M;
LalgITime := Cputime();
SetNthreads(THREADS_TO_BE_USED_IN_THE_LINEAR_ALGEBRA);
sol_Lalg := ModularSolution(M, r ); //: Lanczos:=true);
LalgFTime := Cputime();
Ltime := LalgFTime - LalgITime;
printf "we have found one vector.\n> %o s\n", Ltime;
printf "Does the vector belong to the kernel?\t %o\n", [v_i mod r : v_i in Eltseq(sol_Lalg * Transpose(M))] eq [0 : i in [1 .. N]];

printf "> Saving the vector solution... See file ./outputs/tmp/VECTOR_IN_THE_KERNEL\n";

System("rm -f ./outputs/tmp/VECTOR_IN_THE_KERNEL");
System("touch ./outputs/tmp/VECTOR_IN_THE_KERNEL");
VECTOR_SOLUTION_i := Open("./outputs/tmp/VECTOR_IN_THE_KERNEL", "r+");
pos_i := 0;
for var_i in [1 .. N] do
	Flush(VECTOR_SOLUTION_i);
	Seek(VECTOR_SOLUTION_i, pos_i, 0);
	Put(VECTOR_SOLUTION_i, IntegerToString(sol_Lalg[var_i]) );
	
	pos_i +:= #IntegerToString(sol_Lalg[var_i]);
	Seek(VECTOR_SOLUTION_i, pos_i, 0);
	
	Put(VECTOR_SOLUTION_i, "\n");
	pos_i +:= 1;
end for;

printf "> Computing the dLog!\n";
den := 0;
enu := 0;
for var_i in [1 .. N] do
	den +:= (BTA_lst[var_i,1]*sol_Lalg[var_i]) mod r;
	enu +:= (ALP_lst[var_i,1]*sol_Lalg[var_i]) mod r;
end for;

den := IntegerRing(r)!den;
enu := IntegerRing(r)!enu;
if(den ne 0) then
	LoG_D := -enu/den;
	printf "\nThe DLP was successfully solved it.\n";
	printf "LoG_D * Pt eq Pt_prime? %o\n", IntegerRing()!LoG_D * Pt eq Pt_prime;
	Write("./outputs/tmp/dLog", LoG_D);
else
	printf "\nFailed! Build another associated matrix.\n";
end if;

exit;
