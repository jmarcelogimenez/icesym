%%Script de corrida particionada
%% Datos Boscarol de avance ignición.
%% dtheta proporcional a rpm (como tiempo de combustión constante)
%% Proporción de combustible lineal (a comparar con los valores de Boscarol).

  B = 0.085; #m
  a = 0.044; #m
  l = 0.14;	#m
  Vc = 0.000058; #m^3
  pi = 3.14159;
  Vd=pi*B^2/4*2*a;

CdPlen = 1.0; %1.0; %	%Coeficiente de Descarga del plenum
system(["sed -i 's/CdPlen =.*/CdPlen = " num2str(CdPlen) ";/g' BoscaMono.py"]);

%%Cambio de perfil de leva intake
LvIdat = '"Lv_I-mb.dat"' %%Poner nombre, este es el orig.
system(["sed -i 's/LvI =.*/LvI = loadtxt(" LvIdat ");/g' BoscaMono.py"]);

rpms = [1000,6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7500, 7600,  7700, 7800, 7900, 8000, 8100, 8200];

theta0 = [5.8,5.6985, 5.69, 5.683, 5.67, 5.66, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65]; %%Extrado mapa Boscarol

dtheta = 0.002*2*pi*rpms/60

PHI = 1.25+0.2/2200*[rpms-6000]	%Modelo lineal PHI(rpms)

%%Data Boscarol
powData = [0,215.4, 225.2, 239, 249.4, 255.5, 260.9, 266.9, 272.4, 275.8, 279.4, 282.4, 283.6, 282.6, 283.4, 282.2, 281.3];

torData = [0,25.7, 26, 26.7, 27.1, 26.9, 26.7, 26.6, 26.4, 26.3, 26.3, 26.3, 26, 25.6, 25.4, 24.9, 24.6];
         
for irpm = 1:length(rpms)-2
	system(["sed -i 's/PHI =.*/PHI = " num2str(PHI(irpm)) ";/g' BoscaMono.py"]);
	system(["sed -i 's/rpms =.*/rpms = [" num2str(rpms(irpm)) "];/g' BoscaMono.py"]);
	system(["sed -i 's/theta0 =.*/theta0 = " num2str(theta0(irpm)) ";/g' BoscaMono.py"]);
	system(["sed -i 's/dtheta =.*/dtheta = " num2str(dtheta(irpm)) ";/g' BoscaMono.py"]);
	system(["python monoRun.py >> Mono.log"]); 
end

for irpm = 1:length(rpms)

	  rpm = rpms(irpm); Sp = 4*a*rpm/60;
	  system(["grep '^  0  6' tests/BoscaMono/RPM_" num2str(rpm) "/cyl_0.txt > extractedCyl0.txt"]);    
	  cyl0 = load(["extractedCyl0.txt"]);
	   
	  vol0 = Vc*(1+0.5*8.5*(l/a+1-cos(cyl0(:,3)*pi/180)-((l/a).^2-sin(cyl0(:,3)*pi/180).^2).^(0.5))); %Para 0 y 3

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
end

figure(1),%clf
 hold on, grid
 plot(rpms,power/735.5, "r")
 plot(rpms,powData/4, "r--")
 plot(rpms,torque,"b")
 plot(rpms,torData*9.8/4, "b--")
 xlabel("engine speed [rpm]")
 ylabel("power [CV], torque [Nm]")
 legend({"simulated power", "test power", "simulated torque", "test torque"}, 'location', 'southeast')
 print -dpdf 00-power.pdf 
