%Script para el Post proceso del cilindro

clear
clc

rpms = [ 6000 7000 8000];

  rpm = rpms(1);
  system(["grep '^  0  6' RPM_" num2str(rpm) "/cyl_0.txt > cilindroa.txt"]);    
  cil = load(["cilindroa.txt"]);
 
  rpm = rpms(2);
  system(["grep '^  0  6' RPM_" num2str(rpm) "/cyl_0.txt > cilindrob.txt"]);    
  cilb = load(["cilindrob.txt"]); 
  
  rpm = rpms(3);
  system(["grep '^  0  6' RPM_" num2str(rpm) "/cyl_0.txt > cilindroc.txt"]);    
  cilc = load(["cilindroc.txt"]); 
   
  densidad_max= max(cil(:,5))*1.1;
  presion_max= max(cil(:,6))/1000000*1.1;
  temp_max= max(cil(:,7))*1.1;
  
  x0=[180,180];
  x1=[360,360];
  x2=[540,540];
  y1=[0,densidad_max];
  y2=[0,presion_max];
  y3=[0,temp_max];
   
 figure(1)
 subplot(2,2,1)
 hold on, grid
 plot(cil(:,3), cil(:,5), "b", "LineWidth", 1.5)
 plot(cilb(:,3), cilb(:,5), "g", "LineWidth", 1.5)
 plot(cilc(:,3), cilc(:,5), "c", "LineWidth", 1.5)
 plot(x0 , y1 , "k--")
 plot(x1 , y1 , "k--")
 plot(x2 , y1 , "k--")
 legend("6000","7000", "8000")
 xlabel("Angulo de cigueñal[°]")
 ylabel("Densidad de mezcla [kg/m^3]")
 title ("Densidad de mezcla del cilindro")
 
 subplot(2,2,2)
 hold on, grid
 plot(cil(:,3), cil(:,6)/1000000, "r", "LineWidth", 1.5)
 plot(cilb(:,3), cilb(:,6)/1000000, "g", "LineWidth", 1.5)
 plot(cilc(:,3), cilc(:,6)/1000000, "c", "LineWidth", 1.5)
 plot(x0 , y2 , "k--")
 plot(x1 , y2 , "k--")
 plot(x2 , y2 , "k--")
 legend("6000","7000", "8000")
 xlabel("Angulo de cigueñal[°]")
 ylabel("Presion [MPa]")
 title ("Presion en el cilindro")
 
 subplot (2,2,3)
 hold on, grid
 plot(cil(:,3), cil(:,7), "g", "LineWidth", 1.5)
 plot(cilb(:,3), cilb(:,7), "g", "LineWidth", 1.5)
 plot(cilc(:,3), cilc(:,7), "c", "LineWidth", 1.5)
 plot(x0 , y3 , "k--")
 plot(x1 , y3 , "k--")
 plot(x2 , y3 , "k--")
 xlabel("Angulo de cigueñal[°]")
 ylabel("Temperatura [°K]")
 title ("Temperatura en el cilindro")
 print -dpdf cilindro.pdf
 
 