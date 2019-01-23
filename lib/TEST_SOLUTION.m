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
GAMMA_v := GetingSmoothnessCoefFromStringFile("./outputs/tmp/VECTOR_IN_THE_KERNEL", N);
GAMMA_v := Transpose(Matrix(GAMMA_v));
GAMMA_r := GetingSmoothnessCoefFromStringFile("./outputs/tmp/dLog", 1);

F_time := Cputime();

printf "Does the vector belong to the kernel?\t %o\n", [v_i mod r : v_i in Eltseq(GAMMA_v * Transpose(M))] eq [0 : i in [1 .. N]];
printf "\nThe DLP was successfully solved it.\n";
printf "LoG_D * Pt eq Pt_prime? %o\n", GAMMA_r[1,1] * Pt eq Pt_prime;

exit;
