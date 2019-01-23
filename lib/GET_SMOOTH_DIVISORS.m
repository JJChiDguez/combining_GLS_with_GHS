// ---
// Now, we need to compute (1/LAMBDA)^j for each j = 0, ..., n-1
_, INVERSE_OF_LAMBDA, _ := XGCD(LAMBDA, r);	// (1 / LAMBDA) mod r
POWERS_OF_LAMBDA := [1];
for local_i in [2 .. n ] do
	POWERS_OF_LAMBDA[local_i] := (POWERS_OF_LAMBDA[local_i - 1] * INVERSE_OF_LAMBDA) mod r;
end for;

// ---
PSI_STAR := function(INPUT_DIVISOR)
	PSI_STAR_X := DELTA_1^Degree(INPUT_DIVISOR[1]) * 
				  Evaluate(P_q! [COEF^N : COEF in  ElementToSequence(INPUT_DIVISOR[1])], DELTA_TILDE_1 * w);
	PSI_STAR_Y_PART_u := DELTA_3 * Evaluate(P_q! [COEF^N : COEF in  ElementToSequence(INPUT_DIVISOR[2])], DELTA_TILDE_1 * w);
	PSI_STAR_Y_PART_h := DELTA_4 * Evaluate(P_q! [COEF^N : COEF in  ElementToSequence(h_q mod INPUT_DIVISOR[1])], DELTA_TILDE_1 * w);
	return J_q![PSI_STAR_X, PSI_STAR_Y_PART_u + PSI_STAR_Y_PART_h];
end function;

/* ------------------------------------------------------------------ */
// We clean possible previous data
System("rm -f " cat var_file_value);
System("rm -f " cat var_file_FB_lc);
System("rm -f " cat var_file_R0_lc);
System("rm -f " cat var_file_Alp);
System("rm -f " cat var_file_Bta);

var_R0_cnt := 0;
var_Sp_cnt := 0;
printf "SMOOTH DIVISORS TO BE REACHED: %o\n", var_upbnd;
Log10_bound := Ceiling(Log(10,var_upbnd));
while( cnt0 le var_upbnd ) do
	/* *************************************************************** *
	 * Finding an s-smooth divisor R0                                  *
	 * *************************************************************** */
	var_j := Random(1, var_r);
	R0 := R0 + T[var_j];
	Alp0 := (Alp0 + alpha[var_j]) mod r ;
	Bta0 := (Bta0 + beta[var_j]) mod r ;

	printf "> ";
	while( smooth_test(R0 , smth_bnd) eq false ) do
		var_j := Random(1, var_r);
		R0 := R0 + T[var_j];
		Alp0 := (Alp0 + alpha[var_j]) mod r ;
		Bta0 := (Bta0 + beta[var_j]) mod r ;
	end while;

	/* *************************************************************** *
	 * Updating the factor base Sp, i.e., inserting (without repetion) *
	 * the prime factors of R0 to Sp                                   *
	 * *************************************************************** */	
	
	string_COEFFICIENT := "";
	string_FB_ELEMENTS := "";
	string_R0_ELEMENTS := "";
	string_ALP := Sprintf("%o", Alp0) cat "\n";;
	string_BTA := Sprintf("%o", Bta0) cat "\n";;
			
	D_xs := Factorization(R0[1]);
	SPLIT_OF_R0 := [J_q![D_x[1], R0[2] mod D_x[1]] : D_x in D_xs];			// D_i : irreducible divisor of R0
	SPLIT_COEFF := [D_x[2] : D_x in D_xs];									// c_i : coefficient of each D_i
	SPLIT_SIZE  := #D_xs;													// #(irreducibles divisors of R0)
	
	//SUM_OF_POINTS := J_q!0;
	for i in [1 .. SPLIT_SIZE] do
		ORBIT_UNDER_PSI_STAR := [];												// Sequence that correspond with orbit of Di under PSI_STAR
		ORBIT_UNDER_PSI_STAR[1] := SPLIT_OF_R0[i];								// Di
		for j in [2 .. n] do
			ORBIT_UNDER_PSI_STAR[j] := PSI_STAR(ORBIT_UNDER_PSI_STAR[j - 1]);
		end for;
		
		// Does   Di_j =  PSI_STAR^j(Di)  belong to Sp?
		DOES_IT_BELONG_TIMES_1 := [ Index(Sp, <Di_j[1], Di_j[2]>) : Di_j in ORBIT_UNDER_PSI_STAR ];
		// Does  -Di_j = -PSI_STAR^j(Di)  belong to Sp?
		DOES_IT_BELONG_MINUS_1 := [ Index(Sp, <Di_j[1], (Di_j[2] + h_q) mod Di_j[1]>) : Di_j in ORBIT_UNDER_PSI_STAR ];
				
		ZERO_SEQUENCE := ZeroSequence(IntegerRing(), n);
		if( (DOES_IT_BELONG_TIMES_1 eq ZERO_SEQUENCE) and (DOES_IT_BELONG_MINUS_1 eq ZERO_SEQUENCE) ) then
			var_Sp_cnt +:= 1; indxs := var_Sp_cnt;
			
			Sp[indxs], LOCAL_POSITION := Max([ <Di_j[1], Di_j[2]> : Di_j in ORBIT_UNDER_PSI_STAR ] cat
											         [ <Di_j[1], (Di_j[2] + h_q) mod Di_j[1]> : Di_j in ORBIT_UNDER_PSI_STAR ]);
			
			if(LOCAL_POSITION le n) then
				var_cff :=  SPLIT_COEFF[i] * POWERS_OF_LAMBDA[LOCAL_POSITION];
			else
				LOCAL_POSITION := (LOCAL_POSITION - n);
				var_cff := -SPLIT_COEFF[i] * POWERS_OF_LAMBDA[LOCAL_POSITION] ;
			end if;
			
			Di_lc_string  := IrreducibleDivisorToString([Sp[indxs]], 1);
			string_FB_ELEMENTS := string_FB_ELEMENTS cat Di_lc_string cat "\n";
		else
			if(DOES_IT_BELONG_TIMES_1 ne ZERO_SEQUENCE) then
				indxs, LOCAL_POSITION := Max(DOES_IT_BELONG_TIMES_1);
				var_cff :=  SPLIT_COEFF[i] * POWERS_OF_LAMBDA[LOCAL_POSITION];
			else
				indxs, LOCAL_POSITION := Max(DOES_IT_BELONG_MINUS_1);
				var_cff := -SPLIT_COEFF[i] * POWERS_OF_LAMBDA[LOCAL_POSITION];				
			end if;
		end if;

		//SUM_OF_POINTS +:= J_q![Sp[indxs][1],Sp[indxs][2]] * var_cff;
		
		     Di_lc_string   := IrreducibleDivisorToString([Sp[indxs]], 1);
		string_R0_ELEMENTS  := string_R0_ELEMENTS cat Di_lc_string cat "\n";
		string_COEFFICIENT := string_COEFFICIENT cat Sprintf("%o ", var_cff);
		var_R0_cnt +:= 1;
	end for;
	
	string_COEFFICIENT := string_COEFFICIENT cat "\n";
	
	System("echo -n \"" cat string_COEFFICIENT cat "\" >>" cat var_file_value);
	System("echo -n \"" cat string_FB_ELEMENTS cat "\" >>" cat var_file_FB_lc);
	System("echo -n \"" cat string_R0_ELEMENTS cat "\" >>" cat var_file_R0_lc);
	System("echo -n \"" cat string_ALP cat "\" >>" cat var_file_Alp);
	System("echo -n \"" cat string_BTA cat "\" >>" cat var_file_Bta);
	
	cnt0 +:= 1;
	
	/* ------------------------------------------------ */
	//printf "IS THE SPLIT CORRECT?\t %o\n", SUM_OF_POINTS eq (Alp0*D + Bta0*D_prime);
	mid := Sprintf("%%0%oo/(%%0%oo + 1)", Log10_bound, Log10_bound);
	beg := Sprintf("[%%0%oo]", Log10_bound);
	printf beg cat " #(SMOOTH DIVISORS REACHED)/#(DYNAMIC FACTOR BASE) : " cat mid cat " ~ %2.2o\n", cnt0, cnt0, #Sp, cnt0 / (#Sp + 1)* 1.0;
end while;

Write(var_file_sizes, Sprintf("%o %o %o", cnt0, var_Sp_cnt, var_R0_cnt));
