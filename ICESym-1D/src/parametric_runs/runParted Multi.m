%%Script de corrida particionada
%% Datos Boscarol de avance ignición.
%% dtheta proporcional a rpm (como tiempo de combustión constante)
%% Proporción de combustible lineal (a comparar con los valores de Boscarol).
clear
clc

  B = 0.0865; #m
  a = 0.0445; #m
  l = 0.14;	#m
  Vc = 0.000058; #m^3
  pi = 3.14159;
  Vd=pi*B^2/4*2*a;

rpms = [4000, 5000, 6000, 6400, 6600, 6800, 7000, 7200, 7400, 7500, 7600,  7700, 7800, 7900, 8000, 8100, 8200];

%dt=0; %Variaciones del avance de encendido
%Theta= [5.654866+dt, 5.654866+dt, 5.654866+dt, 5.654866+dt, 5.654866+dt, 5.654866+dt, 5.654866+dt, 5.654866+dt, 5.654866+dt, 5.654866+dt, 5.654866+dt, 5.654866+dt, 5.654866+dt, 5.65+dt, 5.654866+dt, 5.654866+dt, 5.654866+dt]; %%Extraido mapa Boscarol

%dc=1.07453;
%Dtheta = [dc,dc,dc,dc,dc,dc,dc,dc,dc,dc,dc,dc,dc,dc,dc,dc,dc];

%PHI = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]; %1.25+0.2/2200*[rpms-6000]	%Modelo lineal PHI(rpms)

%%Data Boscarol
%powData = [215.4, 225.2, 239, 249.4, 255.5, 260.9, 266.9, 272.4, 275.8, 279.4, 282.4, 283.6, 282.6, 283.4, 282.2, 281.3];

%torData = [25.7, 26, 26.7, 27.1, 26.9, 26.7, 26.6, 26.4, 26.3, 26.3, 26.3, 26, 25.6, 25.4, 24.9, 24.6];

for irpm = 1:length(rpms)
	system(["sed -i 's/rpms =.*/rpms = [" num2str(rpms(irpm)) "];/g'  multicyl4TSI.py"]);
	%system(["sed -i 's/PHI =.*/PHI = " num2str(PHI(irpm)) ";/g'  multicyl4TSI.py"]);
	%system(["sed -i 's/Theta =.*/Theta = " num2str(Theta(irpm)) ";/g'  multicyl4TSI.py"]);
	%system(["sed -i 's/Dtheta =.*/Dtheta = " num2str(Dtheta(irpm)) ";/g'  multicyl4TSI.py"]);
	system(["python multiRun.py >> Multi.log"]); 
end

for irpm = 1:length(rpms)
  rpm = rpms(irpm); Sp = 4*a*rpm/60;
  system(["grep '^  0  5' RPM_" num2str(rpm) "/cyl_0.txt > extractedCyl0.txt"]);    
  cyl0 = load(["extractedCyl0.txt"]);
  
  system(["grep '^  0  5' RPM_" num2str(rpm) "/cyl_1.txt > extractedCyl1.txt"]);    
  cyl1 = load(["extractedCyl1.txt"]);
  
  system(["grep '^  0  5' RPM_" num2str(rpm) "/cyl_2.txt > extractedCyl2.txt"]);    
  cyl2 = load(["extractedCyl2.txt"]);
  
  system(["grep '^  0  5' RPM_" num2str(rpm) "/cyl_3.txt > extractedCyl3.txt"]);    
  cyl3 = load(["extractedCyl3.txt"]);
   
  vol0 = Vc*(1+0.5*8.5*(l/a+1-cos(cyl0(:,3)*pi/180)-((l/a).^2-sin(cyl0(:,3)*pi/180).^2).^(0.5))); %Para 0 y 3
  vol1 = Vc*(1+0.5*8.5*(l/a+1-cos((cyl1(:,3)+540)*pi/180)-((l/a).^2-sin((cyl1(:,3)+540)*pi/180).^2).^(0.5)));
  vol2 = Vc*(1+0.5*8.5*(l/a+1-cos((cyl2(:,3)+180)*pi/180)-((l/a).^2-sin((cyl2(:,3)+180)*pi/180).^2).^(0.5)));
  vol3 = Vc*(1+0.5*8.5*(l/a+1-cos((cyl3(:,3)+360)*pi/180)-((l/a).^2-sin((cyl3(:,3)+360)*pi/180).^2).^(0.5)));
  
  cycle_work = 0;

  cycle_work += trapz(vol0,cyl0(:,6))+trapz(vol1,cyl1(:,6))+trapz(vol2,cyl2(:,6))+trapz(vol3,cyl3(:,6));
  %cycle_work += trapz(vol2,cyl2(:,6))
  totalPower = cycle_work*rpm/60/2;
  totalTorque = cycle_work/(2*pi*2);
  
  pow(irpm) = totalPower/735.5;
  torque(irpm) = totalTorque;
end

  max_pres0=max(cyl0(:,6));
  max_pres1=max(cyl1(:,6));
  max_pres2=max(cyl2(:,6));
  max_pres3=max(cyl3(:,6));

figure(3),%clf
 hold on, grid
 plot(rpms,pow, "b")
 %plot(rpms,torque, "b")
 xlabel("engine speed [rpm]")
 ylabel("power [CV], torque [Nm]")
 legend({"simulated power","simulated torque"}, 'location', 'southeast')
 print -dpdf pow.pdf