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
System("wc -l < ./outputs/tmp/ALP > ./outputs/tmp/N_size");
N := StringToInteger(Gets(Open("./outputs/tmp/N_size", "r")));
M := SparseMatrix(N, N);

I_time := Cputime();
GetingSmoothnessCoefFromStringFile := function(var_File, var_cnt)
	var_ReadFromFile := Open(var_File, "r");
	var_i:= 0; var_k := 1;
	var_SOL :=[];
	for var_l in [1 .. var_cnt] do
		var_sol_aux := Gets(var_ReadFromFile);
		var_sol_aux := StringToIntegerSequence(var_sol_aux);
		var_SOL[var_l] := var_sol_aux;
	end for;

	return var_SOL;
end function;

value_lst := GetingSmoothnessCoefFromStringFile("./outputs/tmp/COEFS", N);
indxs_lst := GetingSmoothnessCoefFromStringFile("./outputs/tmp/INDXS", N);

/* ----------------------------------------------------------------------------- */
printf "[STEP #1] WE READ THE MATRIX FROM THE FILES ./outputs/tmp/COEFS AND ./outputs/tmp/INDXS\n";
time for var_j in [1 .. N] do
	for var_k in [1 .. #value_lst[var_j]] do
		M[var_j, indxs_lst[var_j, var_k]] +:= value_lst[var_j, var_k];
	end for;
end for;

/* ----------------------------------------------------------------------------- */
printf "[STEP #2] WE CREATE TWO STRUCTURES THAT DETERMINE THE TRANPOSE OF THE MATRIX\n";
printf "WE DO IT BY USING THE FILES ./outputs/tmp/COEFS AND ./outputs/tmp/INDXS\n";
I_tr := [];
C_tr := [];
time for var_j in [1 .. N] do
	for var_k in [1 .. #value_lst[var_j]] do
		if( IsDefined(I_tr, indxs_lst[var_j, var_k]) eq false) then
			I_tr[indxs_lst[var_j, var_k]] := [var_j];
			C_tr[indxs_lst[var_j, var_k]] := [value_lst[var_j, var_k]];
		else
			Append(~I_tr[indxs_lst[var_j, var_k]], var_j);
			Append(~C_tr[indxs_lst[var_j, var_k]], value_lst[var_j, var_k]);
		end if;
	end for;
end for;

/* ----------------------------------------------------------------------------- */
printf "[STEP #3] WE CONSTRUCT THE TRANSPOSE OF THE MATRIX BY USING THOSE TWO STRUCTURES\n";
M_TR := SparseMatrix(N, N);
time for var_j in [1 .. #I_tr] do
	if( IsDefined(I_tr, var_j) ) then
		for var_k in [1 .. #I_tr[var_j]] do
			M_TR[var_j, I_tr[var_j, var_k]] +:= C_tr[var_j, var_k];
		end for;
	end if;
end for;

/* ----------------------------------------------------------------------------- */
printf "[STEP #4] WE TRANSPOSE THE MATRIX BY USING THE Transpose() FUNCTION\n";
time tr_M := Transpose(M);
printf "Our transpose matrix constructed is equal to magma transpose output? %o\n", M_TR eq tr_M;

/* ----------------------------------------------------------------------------- */
printf "[STEP #5] WE CREATE THE FILES ./outputs/tmp/tr_INDXS AND ./outputs/tmp/tr_COEFS\n";

F_i := Open("./outputs/tmp/tr_INDXS", "r+");
F_c := Open("./outputs/tmp/tr_COEFS", "r+");

pos_i := 0;
pos_c := 0;
time for var_j in [1 .. #I_tr] do
	if( IsDefined(I_tr, var_j) ) then
		for var_k in [1 .. #I_tr[var_j]] do
			Flush(F_i);
			Flush(F_c);
			//M_TR[var_j, I_tr[var_j, var_k]] := C_tr[var_j, var_k];
			Seek(F_i, pos_i, 0);
			Seek(F_c, pos_c, 0);
			
			Put(F_i, IntegerToString(I_tr[var_j, var_k]) );
			Put(F_c, IntegerToString(C_tr[var_j, var_k]) );
			
			pos_i +:= #IntegerToString(I_tr[var_j, var_k]);
			pos_c +:= #IntegerToString(C_tr[var_j, var_k]);
			
			Seek(F_i, pos_i, 0);
			Seek(F_c, pos_c, 0);

			Put(F_i, " ");
			Put(F_c, " ");
			
			pos_i +:= 1;
			pos_c +:= 1;
		end for;
		
		Seek(F_i, pos_i, 0);
		Seek(F_c, pos_c, 0);	
		
		Put(F_i, "\n");
		Put(F_c, "\n");
		
		pos_i +:= 1;
		pos_c +:= 1;
		
		Flush(F_i);
		Flush(F_c);
	end if;
end for;

/* ----------------------------------------------------------------------------- */
printf "[STEP #6] WE READ THE TRANSPOSE OF THE MATRIX FROM THE FILES ./outputs/tmp/tr_COEFS AND ./outputs/tmp/tr_INDXS\n";
L_TR := SparseMatrix(N, N);
P_TR := SparseMatrix(N, N);

System("wc -l < ./outputs/tmp/tr_INDXS > ./outputs/tmp/L_size");
L := StringToInteger(Gets(Open("./outputs/tmp/L_size", "r")));

vlstr_lst := GetingSmoothnessCoefFromStringFile("./outputs/tmp/tr_COEFS", L);
idxtr_lst := GetingSmoothnessCoefFromStringFile("./outputs/tmp/tr_INDXS", L);

var_cnt := 0;
time for var_j in [1 .. #I_tr] do
	if( IsDefined(I_tr, var_j) ) then
		var_cnt +:= 1;
		for var_k in [1 .. #vlstr_lst[var_cnt]] do
			L_TR[var_cnt, idxtr_lst[var_cnt, var_k]] +:= vlstr_lst[var_cnt, var_k];
		end for;
		
		for var_k in [1 .. #I_tr[var_j]] do
			P_TR[var_j, I_tr[var_j, var_k]] +:= C_tr[var_j, var_k];
		end for;
	end if;
end for;

// L_TR is our transpose matrix loaded from the files tr_COEFS and tr_INDXS
// P_TR should be equal to L_TR;
time L_TR := P_TR;
printf "THE TRANSPOSE OF THE MATRIX CONSTRUCTED BY THE FILES ./outputs/tmp/tr_COEFS AND ./outputs/tmp/tr_INDXS ";
printf "IS EQUAL TO THE MATRIX CONSTRUCTED BY THE FILES ./outputs/tmp/INDXS AND ./outputs/tmp/COEFS? %o\n", Transpose(L_TR) eq M;

exit;
