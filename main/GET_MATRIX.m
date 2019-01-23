/* ------------------------------------------------ */
Nthreads := StringToInteger(Read("./inc/NUMBER_OF_THREADS"));
printf "Creating the %o magma scripts!\n", Nthreads;
var_LINE := "";
for var_i in [1 .. Nthreads] do
	var_LINE := var_LINE cat Sprintf("magma ./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m >> ./outputs/logs/GET_MATRIX_POSITION_USING_THREAD_%o.log & ", var_i, var_i);
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), "clear;" : Overwrite:=true);
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), "cnt0 := 0;");
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), Sprintf("var_FB_size := \"./outputs/logs/var_FB_size_%o\";", var_i));
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), "load \"./inc/FACTOR_BASE.m\";");
	
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), "print \"\\nThe main loop will start \";");
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), "printf \"\\nObtaining positions in the matrix... \";");

	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), Sprintf("var_file_value := \"./outputs/tmp/values_%o\";", var_i));
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), Sprintf("var_file_R0_lc := \"./outputs/tmp/R0_lc_%o\";", var_i));
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), Sprintf("var_file_sizes := \"./outputs/tmp/sizes_%o\";", var_i));
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), Sprintf("var_file_indxs := \"./outputs/tmp/indxs_%o\";", var_i));
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), Sprintf("var_file_value := \"./outputs/tmp/values_%o\";", var_i));
	
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), Sprintf("var_file_ALP := \"./outputs/tmp/Alp_cff_%o\";", var_i));
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), Sprintf("var_file_BTA := \"./outputs/tmp/Bta_cff_%o\";", var_i));
	
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), "RltsITime := Cputime();");
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), "load \"./lib/COEFFICIENTS.m\";");
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), "RltsFTime := Cputime();");
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), "Rtime := RltsFTime - RltsITime;");
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), "printf \"... all the positions were saved.\\n\";");
	
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), Sprintf("Write(\"./outputs/logs/Timing_R0s_%o\", ", var_i) cat "Sprintf(\"%o\", Rtime));");
	Write(Sprintf("./outputs/GET_MATRIX_POSITION_USING_THREAD_%o.m", var_i), "exit;");
end for;

main_string := "echo \"" cat var_LINE cat " \" > ./SCRIPT_GET_MATRIX_POSITION.sh";
System(main_string);
System("chmod +x ./SCRIPT_GET_MATRIX_POSITION.sh");
exit;
