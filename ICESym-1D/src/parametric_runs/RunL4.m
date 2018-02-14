%%Script de varias corridas parametrizadas por L4
clear
clc

L4tests = [0.15 0.25 0.35]; %Prolongacion de escape, original = 0.25
for iL4 = 1:length(L4tests)
L4 = L4tests(iL4);
cmd = ["sed -i 's/L4 =.*/L4 = " num2str(L4) ";/g' BoscaMono.py"]
system(cmd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rpms = [6000, 6200, 6400,6600, 6800, 7000, 7200, 7400, 7500, 7600,  7700, 7800, 7900, 8000, 8100, 8200];
theta0 = [ 5.6985, 5.69, 5.683, 5.67, 5.66, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65,5.65,5.65];

for irpm = 1:length(rpms)
  system(["sed -i 's/rpms =.*/rpms = [" num2str(rpms(irpm)) "];/g' BoscaMono.py"]);
  system(["sed -i 's/theta0 =.*/theta0 = " num2str(theta0(irpm)) ";/g' BoscaMono.py"]);
  system(["python monoRun.py >> Mono.log"]); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmd = ["mv tests/BoscaMono tests/BoscaMono_L4_" num2str(L4)]
system(cmd);
end
