/* ***************************************************************************** *
 * ----------------------------------------------------------------------------- *
 * Positive integers to be used:                                                 *
 *                                                                               *
 * n > 0 and l > 0 s.t. (n,l) = 1,                                               *
 * q <- 2^n,                                                                     *
 * N <- 2^l.                                                                     *
 * ----------------------------------------------------------------------------- *
 * Binary finite fields to be used:                                              *
 *                                                                               *
 * F_2 will denote the finite field GF(2),                                       *
 * F_q will denote the degree m extension field of F_2,                          *
 * K_q will denote the degree n extension field of F_2,                          *
 * F_qn will denote the degree n extension field of F_q.                         *
 * ----------------------------------------------------------------------------- *
 * Polynomial rings to be used:                                                  *
 *                                                                               *
 * P_2  <- F_2[t],                                                               *
 * P_q  <- F_q[w],                                                               *
 * P_qn <- F_2[z].                                                               *
 * ----------------------------------------------------------------------------- *
 * ***************************************************************************** */

n := 2; l := 8191; q := 2^n; N := 2^l;

F_2 := GF(2);
P_2<t> := PolynomialRing(F_2);

poly_q_n := P_2!IrreduciblePolynomial(GF(2), n);
poly_q_l := P_2!IrreduciblePolynomial(GF(2), l);

F_q<u> := ext<F_2| poly_q_n>;
K_q<s> := ext<F_2| poly_q_l>;
F_qn<v>:= ext<F_q| poly_q_l>;

P_q<w> := PolynomialRing(F_q);
P_qn<z> := PolynomialRing(F_qn);
