%%Script de varias corridas parametrizadas por L4

L4tests = [0.25 0.35 0.5]; %0.35; %	%Prolongacion de escape, original = 0.25
for iL4 = 1:length(L4tests)
L4 = L4tests(iL4);
cmd = ["sed -i 's/L4 =.*/L4 = " num2str(L4) ";/g' BoscaMono.py"]
system(cmd);
runParted
cmd = ["mv tests/BoscaMono tests/BoscaMono_L4_" num2str(L4)]
system(cmd);
end
