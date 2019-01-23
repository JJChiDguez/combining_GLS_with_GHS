clear;

/* ----------------------------------------------------------------------------- */
load "./inc/FINITE_FIELD.m";

/* ----------------------------------------------------------------------------- */
load "./inc/tmp/ELLIPTIC_CURVE_parameters";
E_qn := EllipticCurve([F_qn| 1, a_qn, 0, 0, b_qn]);
E_qn_size := c*r;

// We randomly select two order-r points
Pt := Random(E_qn) * c;
Pt_prime := Random(E_qn) * c;

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
H_q, W_q := WeilDescent(E_qn, F_q, F_qn!1);
f_q, h_q := HyperellipticPolynomials(H_q);
J_q := Jacobian(H_q);

D := W_q(Pt);
D_prime := W_q(Pt_prime);


// ---
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


DELTA_TILDE_1, DELTA_1, DELTA_3, DELTA_4 := GET_deltas_134(D);

// Once the coeffients DELTA_TILDE_1, DELTA_1, DELTA_3, and DELTA_4 are given, then 
// we can define the PSI_STAR on the genus-g hyperelliptic curve and its jacobian. 

// ---
PSI_STAR := function(INPUT_DIVISOR)
	PSI_STAR_X := DELTA_1^Degree(INPUT_DIVISOR[1]) * 
				  Evaluate(P_q! [COEF^N : COEF in  ElementToSequence(INPUT_DIVISOR[1])], DELTA_TILDE_1 * w);
	PSI_STAR_Y_PART_u := DELTA_3 * Evaluate(P_q! [COEF^N : COEF in  ElementToSequence(INPUT_DIVISOR[2])], DELTA_TILDE_1 * w);
	PSI_STAR_Y_PART_h := DELTA_4 * Evaluate(P_q! [COEF^N : COEF in  ElementToSequence(h_q mod INPUT_DIVISOR[1])], DELTA_TILDE_1 * w);
	return J_q![PSI_STAR_X, PSI_STAR_Y_PART_u + PSI_STAR_Y_PART_h];
end function;

STRING_0 := Sprintf("%%0%oo", Ceiling(Log(10, r)));

printf "Now, we test the endomorphism PSI_STAR!\n";
for i in [1 .. 100] do
	R := Random(2, r-1); 
	BOOL := (LAMBDA*R*D) eq PSI_STAR(R*D);
	if(BOOL eq false) then printf "[ERROR]\n"; exit; end if;
	printf "(%03o) PSI_STAR((" cat STRING_0 cat ") * D) eq (" cat STRING_0 cat " * " cat STRING_0 cat ") * D? %o,\n", i, R, LAMBDA, R, BOOL;
end for;

i := 1; Dt := PSI_STAR(D);
while(Dt ne D) do
	i +:= 1; Dt := PSI_STAR(Dt);
end while;
printf "\n(PSI_STAR)^ %o = 1\n", i;

//exit;
