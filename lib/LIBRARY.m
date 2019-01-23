/* ***************************************************************************** *
 * ***************************************************************************** */


/* Naive smoothness test implementation                                          */
smooth_test := function(D, smth_bnd)
	var_aux := Factorisation (D[1]);
	var_AuX := Degree(var_aux[#var_aux, 1]);
	if(var_AuX gt smth_bnd) then
		return false;
	end if;

	return true;
end function;
 

/* Saving s-smooth divisors in a file                                            */
IrreducibleDivisorToString := function(var_D_prime, var_cnt)
	var_SOL := "";
	for var_k in [1 .. var_cnt] do
		var_u := ElementToSequence(var_D_prime[var_k, 1]);
		var_len_u := #var_u;
		var_v := ElementToSequence(var_D_prime[var_k, 2]);
		var_len_v := #var_v;
		line_len := Sprintf("%o %o\n", var_len_u, var_len_v);
		String_u := "";
		for var_J in [1 .. var_len_u] do
			var_ui := ElementToSequence(var_u[var_J]);
			var_len_ui := #var_ui;
			line_ui := "";
			for var_K in [1 .. (var_len_ui - 1)] do
				line_ui := line_ui cat Sprintf("%o ", var_ui[var_K]);
			end for;
			line_ui := line_ui cat Sprintf("%o", var_ui[var_len_ui]);
			String_u := String_u cat line_ui cat "\n";
		end for;
		String_v := "";
		for var_J in [1 .. var_len_v] do
			var_vi := ElementToSequence(var_v[var_J]);
			var_len_vi := #var_vi;
			line_vi := "";
			for var_K in [1 .. (var_len_vi - 1)] do
				line_vi := line_vi cat Sprintf("%o ", var_vi[var_K]);
			end for;
			line_vi := line_vi cat Sprintf("%o", var_vi[var_len_vi]);
			String_v := String_v cat line_vi cat "\n";
		end for;
		var_SOL := var_SOL cat (line_len cat String_u cat String_v);
	end for;
	return var_SOL;
end function;


/* Saving s-smooth divisors with its coefficients in a file                     */
GoodRelationWithCoefsToString := function(var_D_prime, var_Alp, var_Bta, var_cnt)
	var_SOL := "";
	for var_k in [1 .. var_cnt] do
		var_u := ElementToSequence(var_D_prime[var_k, 1]);
		var_len_u := #var_u;
		var_v := ElementToSequence(var_D_prime[var_k, 2]);
		var_len_v := #var_v;
		line_len := Sprintf("%o %o\n", var_len_u, var_len_v);
		line_Alp:= Sprintf("%o\n", var_Alp[var_k]);
		line_Bta:= Sprintf("%o\n", var_Bta[var_k]);
		String_u := "";
		for var_J in [1 .. var_len_u] do
			var_ui := ElementToSequence(var_u[var_J]);
			var_len_ui := #var_ui;
			line_ui := "";
			for var_K in [1 .. (var_len_ui - 1)] do
				line_ui := line_ui cat Sprintf("%o ", var_ui[var_K]);
			end for;
			line_ui := line_ui cat Sprintf("%o", var_ui[var_len_ui]);
			String_u := String_u cat line_ui cat "\n";
		end for;
		String_v := "";
		for var_J in [1 .. var_len_v] do
			var_vi := ElementToSequence(var_v[var_J]);
			var_len_vi := #var_vi;
			line_vi := "";
			for var_K in [1 .. (var_len_vi - 1)] do
				line_vi := line_vi cat Sprintf("%o ", var_vi[var_K]);
			end for;
			line_vi := line_vi cat Sprintf("%o", var_vi[var_len_vi]);
			String_v := String_v cat line_vi cat "\n";
		end for;
		var_SOL := var_SOL cat (line_Alp cat line_Bta cat line_len cat String_u cat String_v);
	end for;
	return var_SOL;
end function;


 
/* Geting divisors from File                                                     */
GetingIrreducibleDivisorsFromStringFile := function(var_File, var_cnt)
	var_ReadFromFile := Open(var_File, "r");
	var_i:= 0; var_k := 1;
	var_SOL :=[];
	for var_l in [1 .. var_cnt] do
		var_uv_len := Gets(var_ReadFromFile);
		while(var_uv_len eq "") do var_uv_len := Gets(var_ReadFromFile); end while;
		var_uv_len := StringToIntegerSequence(var_uv_len);
		var_u_len := var_uv_len[1];
		var_v_len := var_uv_len[2];
		var_aux_ui := [];
		for var_J in [1 .. var_u_len] do
			Append(~var_aux_ui, F_q!StringToIntegerSequence(Gets(var_ReadFromFile)));
		end for;
	
		var_aux_vi := [];
		for var_J in [1 .. var_v_len] do
			Append(~var_aux_vi, F_q!StringToIntegerSequence(Gets(var_ReadFromFile)));
		end for;
		
		Append(~var_SOL, J_q![P_q!var_aux_ui, P_q!var_aux_vi]);
	end for;

	return var_SOL;
end function;


/* Geting divisors with coefficients from File                                   */
GetingGoodRelationWithCoefFromStringFile := function(var_File, var_cnt)
	var_ReadFromFile := Open(var_File, "r");
	var_i:= 0; var_k := 1;
	var_SOL :=[]; var_Alp := []; var_Bta := [];
	for var_l in [1 .. var_cnt] do
		var_Alp_i := Gets(var_ReadFromFile);
		while(var_Alp_i eq "") do var_Alp_i := Gets(var_ReadFromFile); end while;
		var_Alp_i := StringToInteger(var_Alp_i);
		
		var_Bta_i := StringToInteger(Gets(var_ReadFromFile));
		var_uv_len := StringToIntegerSequence(Gets(var_ReadFromFile));
		var_u_len := var_uv_len[1];
		var_v_len := var_uv_len[2];
		var_aux_ui := [];
		for var_J in [1 .. var_u_len] do
			Append(~var_aux_ui, F_q!StringToIntegerSequence(Gets(var_ReadFromFile)));
		end for;
	
		var_aux_vi := [];
		for var_J in [1 .. var_v_len] do
			Append(~var_aux_vi, F_q!StringToIntegerSequence(Gets(var_ReadFromFile)));
		end for;
		
		Append(~var_SOL, J_q![P_q!var_aux_ui, P_q!var_aux_vi]);
		Append(~var_Alp, var_Alp_i);
		Append(~var_Bta, var_Bta_i);
	end for;

	return var_SOL, var_Alp, var_Bta;
end function;


/* Geting coefficients Alp, Bta, M[i,j], or j from file                          */
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

