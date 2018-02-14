%Script rápido para el calculo de potencia y torque

rpms = [6000, 6400,6600, 6800, 7000, 7200, 7400, 7500, 7600, 7700, 7800, 7900, 8000];%

  B = 0.085; #m
  a = 0.044; #m
  l = 0.14;	#m
  Vc = 0.000058; #m^3
  pi = 3.14159;
  Vd=pi*B^2/4*2*a;
         
         
for irpm = 1:length(rpms)

  rpm = rpms(irpm); Sp = 4*a*rpm/60;
  system(["grep '^  0  6' RPM_" num2str(rpm) "/cyl_0.txt > extractedCyl0.txt"]);

  cyl0 = load(["extractedCyl0.txt"]);
  vol0 = Vc*(1+0.5*8.5*(l/a+1-cos(cyl0(:,3)*pi/180)-((l/a).^2-sin(cyl0(:,3)*pi/180).^2).^(0.5)));

  cycle_work = 0;

  cycle_work += trapz(vol0,cyl0(:,6));

  totalPower = cycle_work*rpm/60/2;
  totalTorque = cycle_work/(2*pi*2);
  
  pcyl_max(irpm) = max(cyl0(:,6));
  imep(irpm)= cycle_work/Vd;
  fmep(irpm) = 1e5*(0.4+0.005e-5*pcyl_max(irpm)+0.09*Sp+9e-4*Sp^2); %From Lopez
  mecEff(irpm) = 1 - fmep(irpm)/imep(irpm); 
  
  power(irpm) = totalPower;
  power(irpm) = power(irpm)*mecEff(irpm);
  torque(irpm) = totalTorque;
  torque(irpm) = torque(irpm)*mecEff(irpm);

figure(2)
clf
 hold on, grid
 plot(vol0*1000000,cyl0(:,6)/1000000, "r")
xlabel("Volume[cm^3]")
ylabel("Pressure [MPa]")
legend("PV diagram")
 print -dpdf cycle.pdf
system(["mv cycle.pdf cycle" num2str(rpm) ".pdf"]);
end

system(["rm extractedCyl*.txt"]);
system(["mv *.pdf postProc*"]);

fmep
