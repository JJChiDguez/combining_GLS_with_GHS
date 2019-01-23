clear;

/* ----------------------------------------------------------------------------- */
load "./inc/FINITE_FIELD.m";
/* ----------------------------------------------------------------------------- */
load "./inc/tmp/ELLIPTIC_CURVE_parameters";
E_qn := EllipticCurve([F_qn| 1, a_qn, 0, 0, b_qn]);
E_qn_size := c*r;

load "./inc/tmp/ELLIPTIC_CURVE_instance";
load "./inc/SMOOTHNESS_PARAMETERS.m";

// Frobenius morphisms
pi_2n := function(var_P)
	return [var_P[1]^N, var_P[2]^N];
end function;

// Polynomial ring overe F_r
P_r<T> := PolynomialRing(GF(r));

// Computing the isomorphism from C_qn into E_qn
delta := F_q!0;
if( (n*l mod 2) eq 1) then
	g_cyc := a_qn + a_qn^N;
	for i in [0 .. ( (n*l -1) div 2 )] do
		delta +:= g_cyc^(2^(2*i) mod (q-1));
	end for;
else
	delta := Roots(w^2 + w + F_q!a_qn^N + F_q!a_qn)[1,1];	
end if;

// The isomorphism is given by (x,y) -> (x, y + delta*x)
phi := function(v_Pt, v_delta)
	return [v_Pt[1], v_Pt[2] + v_delta*v_Pt[1]];
end function;

// Now, the GLS endomorphism is given by the composition  phi * (pi^l)
psi := function(v_Pt, delta)
	return E_qn!phi(pi_2n(v_Pt), delta);
end function;
	
// Finding the eingevalue lambda of the GLS endomorphism
print "finding the eingevalue lambda of the endomorphism psi!";
Qt := psi(Pt, delta);
for el_i in  Roots(T^(2*n) - 1) do
	LAMBDA := IntegerRing()!el_i[1];
	LAMBDA_TIMES_Pt := LAMBDA * Pt;
	if( (LAMBDA_TIMES_Pt eq Qt) and (LAMBDA_TIMES_Pt[1] ne Pt[1]) )then
		break;
	end if;
end for;


/* ----------------------------------------------------------------------------- */
load "./inc/tmp/GENUS_g_HCURVE_parameters";
H_q := HyperellipticCurve(f_q, h_q);
J_q := Jacobian(H_q);
load "./inc/tmp/GENUS_g_HCURVE_instance";

// -----------------------------------------
GET_deltas_134 := function(LOCAL_D)

	degree_u:= Degree(LOCAL_D[1]);						// deg u
	u_local := [COEF^N : COEF in Eltseq(LOCAL_D[1])];	// Coefficients of (^\sigma)u
	u_prime := Eltseq((LAMBDA*LOCAL_D)[1]);				// Coefficients of u'

	degree_v:= Degree(LOCAL_D[2]);						//deg v
	v_local := [COEF^N : COEF in Eltseq(LOCAL_D[2])];	// Coefficients of (^\sigma)v
	v_prime := Eltseq((LAMBDA*LOCAL_D)[2]);				// Coefficients of v'	

	h_local := [COEF^N : COEF in Eltseq(h_q mod LOCAL_D[1])];	// Coefficients of (^\sigma)(h mod u)

	// computing DELTA_TILDE_1
	DENOMINATOR := F_q!0;
	repeat
		local_j := Random(1, degree_u);	local_i := local_j + 1;
	until ( (u_local[local_i] * u_prime[local_j]) ne F_q!0 );
	LOCAL_DELTA_TILDE_1 := (u_local[local_j] * u_prime[local_i]) / (u_local[local_i] * u_prime[local_j]);

	LOCAL_DELTA_1 := 1 / LOCAL_DELTA_TILDE_1;						// DELTA_TILDE_1^{-1}
	DELTA_i := [LOCAL_DELTA_1^(i-1) : i in [1 .. (degree_v + 1)] ];	// DELTA_1^i for i=0, ..., degree_v

	// computing DELTA_3 and DELTA_4
	repeat
		local_i := Random([1 .. (degree_v + 1)]);
		local_j := Random([1 .. (degree_v + 1)]);

	until ( (v_local[local_i] * (h_local[local_i] * v_local[local_j] + h_local[local_j] * v_local[local_i])) ne F_q!0 );

	LOCAL_DELTA_4 :=	(v_prime[local_i] * DELTA_i[local_i] * v_local[local_j] + 
						 v_prime[local_j] * DELTA_i[local_j] * v_local[local_i]
						)/
						(h_local[local_i] * v_local[local_j] + h_local[local_j] * v_local[local_i]);

	LOCAL_DELTA_3 := (v_prime[local_i] * DELTA_i[local_i] + LOCAL_DELTA_4 * h_local[local_i]) / v_local[local_i];


	return LOCAL_DELTA_TILDE_1, LOCAL_DELTA_1, LOCAL_DELTA_3, LOCAL_DELTA_4;

end function;

/* ----------------------------------------------------------------------------- */
printf "Computing the coefficients delta_1, delta_3, delta_4 and delta_5\n";
time DELTA_TILDE_1, DELTA_1, DELTA_3, DELTA_4 := GET_deltas_134(D);

// -------------------------
// Once the coeffients DELTA_TILDE_1, DELTA_1, DELTA_3, and DELTA_4 are given, then 
// we can define the PSI_STAR on the genus-g hyperelliptic curve and its jacobian. 

/* ----------------------------------------------------------------------------- */
load "./lib/LIBRARY.m";
load "./lib/RANDOM_WALK.m";

/* ----------------------------------------------------------------------------- */
Nthreads := StringToInteger(Read("./inc/NUMBER_OF_THREADS"));
var_acc := 0;
for var_i in [1 .. smth_bnd] do
	// #[div(u,v) such that u is irreducible and of degree var_i] 
	// ~ #[Irreducible polynomials in F_q] / 2. However, because 
	// we can reduced by the action of PSI_STAR, we have
	// #[div(u,v) such that u is irreducible and of degree var_i]
	// ~ #[Irreducible polynomials in F_q] / (2*n)).
	
	var_acc +:= NumberOfPrimePolynomials(q, var_i) div (2*n);
end for;

// Each thread will find var_acc / Nthreads;
var_upbnd := (var_acc div Nthreads) + 1;
printf "Creating the %o magma scripts!\n", Nthreads;

/* ----------------------------------------------------------------------------- */
var_LINE := "";
for var_i in [1 .. Nthreads] do
	var_LINE := var_LINE cat Sprintf("magma ./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m >> ./outputs/logs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.log & ", var_i, var_i);
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "clear;" : Overwrite:=true);
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "load \"./inc/FRAMEWORK.m\";");
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "cnt0 := 0;");
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), Sprintf("var_upbnd := %o;", var_upbnd));
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), Sprintf("var_r := %o; smth_bnd := %o; LAMBDA := %o;", var_r, smth_bnd, LAMBDA));
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), Sprintf("DELTA_TILDE_1 := %o; DELTA_1 := %o; DELTA_3 := %o; DELTA_4 := %o;", DELTA_TILDE_1, DELTA_1, DELTA_3, DELTA_4));
	
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "print \"\\nThe main loop will start \";");
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "printf \"\\nThe relation phase is started... \";");

	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), Sprintf("var_file_FB_lc := \"./outputs/tmp/FB_lc_%o\";", var_i));	
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), Sprintf("var_file_Alp   := \"./outputs/tmp/Alp_cff_%o\";", var_i));
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), Sprintf("var_file_Bta   := \"./outputs/tmp/Bta_cff_%o\";", var_i));
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), Sprintf("var_file_value := \"./outputs/tmp/values_%o\";", var_i));
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), Sprintf("var_file_R0_lc := \"./outputs/tmp/R0_lc_%o\";", var_i));
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), Sprintf("var_file_sizes := \"./outputs/tmp/sizes_%o\";", var_i));	
	
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "RltsITime := Cputime();");
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "load \"./lib/GET_SMOOTH_DIVISORS.m\";");
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "RltsFTime := Cputime();");
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "Rtime := RltsFTime - RltsITime;");
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), Sprintf("Write(\"./outputs/logs/R0_cnt_%o\", ", var_i) cat "Sprintf(\"%o\", cnt0));");
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), Sprintf("Write(\"./outputs/logs/Sp_cnt_%o\", ", var_i) cat "Sprintf(\"%o\", #Sp));");
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "printf \"... the relation phase is finished.\\n\";");
	
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "printf \"\\nThe dynamic factor base size is %o\", #Sp;");
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "line_s:= Sprintf(\"Rtime := %o;\", Rtime);");
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), Sprintf("Write(\"./outputs/logs/Timing_R0s_%o\", ", var_i) cat "Sprintf(\"%o\", Rtime));");
	Write(Sprintf("./outputs/GET_SMOOTH_DIVISORS_USING_THREAD_%o.m", var_i), "exit;");
end for;

/* ----------------------------------------------------------------------------- */
main_string := "echo \"" cat var_LINE cat " \" > ./SCRIPT_GET_SMOOTH_DIVISORS.sh";
System(main_string);
System("chmod +x ./SCRIPT_GET_SMOOTH_DIVISORS.sh");
exit;
