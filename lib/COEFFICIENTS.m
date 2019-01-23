var_size_AuX := StringToIntegerSequence(Gets(Open(var_file_sizes, "r")));
var_V := GetingSmoothnessCoefFromStringFile(var_file_value, var_size_AuX[1]);

var_file_ALP, var_file_BTA;
var_A := GetingSmoothnessCoefFromStringFile(var_file_ALP, var_size_AuX[1]);
var_B := GetingSmoothnessCoefFromStringFile(var_file_BTA, var_size_AuX[1]);

R0_lc := GetingIrreducibleDivisorsFromStringFile(var_file_R0_lc, var_size_AuX[3]);

main_aux_string := "";
for var_i in [ 1 .. var_size_AuX[1]] do
	var_len_fcts := #var_V[var_i];
	/* ------------------------------------------------ */
	var_string_aux := "";
	printf "ROW #%o of %o\t", var_i, var_size_AuX[1];
	//ACC := J_q!0;
	for var_k in [1 .. var_len_fcts] do
		cnt0 +:= 1;
		var_aux_Indxs := Index(FB, R0_lc[cnt0]);
		var_string_aux := var_string_aux cat Sprintf("%o ", var_aux_Indxs);
		//ACC +:= FB[var_aux_Indxs] * var_V[var_i][var_k];
	end for;
	//printf "IS THE SPLIT CORRECT?\t %o\t", ACC eq (var_A[var_i][1]*D + var_B[var_i][1]*D_prime);
	main_aux_string   := main_aux_string cat Sprintf("%o", var_string_aux) cat "\n";

	if ( ( (var_i mod 27) eq 0 ) or (var_i eq var_size_AuX[1]) ) then
		System("echo -n \"" cat main_aux_string cat "\" >>" cat var_file_indxs);
		main_aux_string := "";
		printf "ROWS saved!\n";
	else 
		printf "\n";
	end if;
	/* ------------------------------------------------ */
end for;
