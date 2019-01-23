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
Sp := []; T := []; alpha := []; beta := [];

T, alpha, beta := GetingGoodRelationWithCoefFromStringFile("./outputs/RANDOM_WALK", var_r);
R0_var, Alp0_var, Bta0_var := GetingGoodRelationWithCoefFromStringFile("./outputs/T0", 1);

Alp0 := Alp0_var[1];
Bta0 := Bta0_var[1];
R0 := R0_var[1];
