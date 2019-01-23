alpha := Random(1,r-1);
beta  := Random(1,r-1);
T     := alpha*D + beta*D_prime;
String_Random_Walk := GoodRelationWithCoefsToString([T], [alpha], [beta], 1);
Write("./outputs/RANDOM_WALK", String_Random_Walk : Overwrite:=true);
	
for i in [2 .. var_r] do
	alpha := Random(1,r-1);
	beta  := Random(1,r-1);
	T     := alpha*D + beta*D_prime;
	String_Random_Walk := GoodRelationWithCoefsToString([T], [alpha], [beta], 1);
	Write("./outputs/RANDOM_WALK", String_Random_Walk);
end for;

Alp0 := Random(1,r-1);
Bta0 := Random(1,r-1);

R0 := Alp0*D + Bta0*D_prime;;
String_T0 := GoodRelationWithCoefsToString([R0], [Alp0], [Bta0], 1);
Write("./outputs/T0", String_T0 : Overwrite:=true);
