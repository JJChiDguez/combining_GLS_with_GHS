a_qn := F_qn!u;
b_qn := v^30 + v^28 + v^27 + v^26 + v^25 + v^22 + v^20 + v^19 + v^18 + v^13 + v^12 + v^10 + v^9 + v^8 + v^7 + v^6 + v^5 + v^3 + 1;

line_s:= Sprintf("a_qn := %o;", a_qn);
Write("./inc/tmp/ELLIPTIC_CURVE_parameters", line_s : Overwrite:=true);
line_s:= Sprintf("b_qn := %o;", b_qn);
Write("./inc/tmp/ELLIPTIC_CURVE_parameters", line_s);

E_qn := EllipticCurve([F_qn| 1, a_qn, 0, 0, b_qn]);

E_qn_size := #E_qn;
size_fcts := Factorization(E_qn_size);

r := size_fcts[#size_fcts,1];
c := E_qn_size div r;

pi_x := v^355 / v^133 + u;
pi_Y := Roots(z^2 + z * pi_x + pi_x^3 + a_qn * (pi_x^2) + b_qn);
pi_y := pi_Y[1,1];

Pt_prime := c * E_qn![pi_x, pi_y];
Pt := c * Random(E_qn);

line_s:= Sprintf("r := %o;", r);
Write("./inc/tmp/ELLIPTIC_CURVE_parameters", line_s);
line_s:= Sprintf("c := %o;", c);
Write("./inc/tmp/ELLIPTIC_CURVE_parameters", line_s);

line_s:= Sprintf("Pt := E_qn!%o;", [Pt[1], Pt[2]]);
Write("./inc/tmp/ELLIPTIC_CURVE_instance", line_s : Overwrite:=true);
line_s:= Sprintf("Pt_prime := E_qn!%o;", [Pt_prime[1], Pt_prime[2]]);
Write("./inc/tmp/ELLIPTIC_CURVE_instance", line_s);
