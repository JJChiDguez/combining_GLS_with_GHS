H_q, W_q := WeilDescent(E_qn, F_q, F_qn!1);
f_q, h_q := HyperellipticPolynomials(H_q);

line_s:= Sprintf("h_q := %o;", P_q!h_q);
Write("./inc/tmp/GENUS_g_HCURVE_parameters", line_s : Overwrite:=true);
line_s:= Sprintf("f_q := %o;", P_q!f_q);
Write("./inc/tmp/GENUS_g_HCURVE_parameters", line_s);

J_q := Jacobian(H_q);

D := W_q(Pt);
D_prime := W_q(Pt_prime);

line_s:= Sprintf("D := J_q!%o;", [D[1], D[2]]);
Write("./inc/tmp/GENUS_g_HCURVE_instance", line_s : Overwrite:=true);
line_s:= Sprintf("D_prime := J_q!%o;", [D_prime[1], D_prime[2]]);
Write("./inc/tmp/GENUS_g_HCURVE_instance", line_s);
