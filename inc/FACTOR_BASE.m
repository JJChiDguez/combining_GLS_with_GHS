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

/* ----------------------------------------------------------------------------- */
FB := Set([]); //IndexedSet([]);
M_ij := [];

printf "\nREADING ALL THE DYNAMIC FACTOR BASE - step 1\n";
col_size := 0; row_size := 0;
Nthreads := StringToInteger(Read("./inc/NUMBER_OF_THREADS"));
I_time := Cputime();
for var_i in [1 .. Nthreads] do
	i_time := Cputime();
	var_name := Sprintf("./outputs/tmp/FB_lc_%o", var_i);
	var_size := StringToIntegerSequence(Gets(Open(Sprintf("./outputs/tmp/sizes_%o", var_i), "r")));
	FB_aux := GetingIrreducibleDivisorsFromStringFile(var_name, var_size[2]);
	col_size +:= #FB_aux;
	row_size +:= var_size[1];
	FB := FB join Set(FB_aux); //IndexedSet(FB_aux);
	f_time := Cputime();
	printf "We have read %o local factor bases, and this iteration takes %os...\n", var_i, f_time - i_time;
end for;
F_time := Cputime();

printf "\n Total time used %o s. Now... We can build the associated matrix\n", F_time - I_time;
FB := SetToSequence(FB);
FB_size:= #FB;
printf "\nThe Dynamic Factor Base size is %o, and  it has %o bits.\n", FB_size, Log(2, FB_size);
Write(var_FB_size, Sprintf("%o", FB_size));
