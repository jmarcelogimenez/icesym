%Script para el Post proceso del tubo 1

clear
clc

rpms = [ 6000 7000 8000];

  rpm = rpms(1);
  system(["grep '^ 52  3' RPM_" num2str(rpm) "/tube_1.txt > tuboadm1a.txt"]);    
  tubo1a = load(["tuboadm1a.txt"]);
  
  rpm = rpms(2);
  system(["grep '^ 52  3' RPM_" num2str(rpm) "/tube_1.txt > tuboadm2b.txt"]);    
  tubo1b = load(["tuboadm2b.txt"]);
  
  rpm = rpms(3);
  system(["grep '^ 52  3' RPM_" num2str(rpm) "/tube_1.txt > tuboadm3c.txt"]);    
  tubo1c = load(["tuboadm3c.txt"]);
   
  velocidad_max= max(tubo1c(:,6));
  velocidad_min= min(tubo1c(:,6));
  densidad_max= max(tubo1c(:,5));
  densidad_min= min(tubo1c(:,5));
  presion_max= max(tubo1c(:,7));
  presion_min= min(tubo1c(:,7));
  
  punto = load(["puestapto.txt"]);
    
  IVO = punto(1)
  IVC = punto(2)
  EVO = punto(3)
  EVC = punto(4)
  
  x0=[IVO/pi*180,IVO/pi*180];
  x1=[IVC/pi*180,IVC/pi*180];
  x2=[EVO/pi*180,EVO/pi*180];
  x3=[EVC/pi*180,EVC/pi*180];
  y1=[velocidad_min,velocidad_max];
  y2=[densidad_min,densidad_max];
  y3=[presion_min,presion_max];  
   
 figure(1),%clf
 hold on, grid
 plot(tubo1a(:,3), tubo1a(:,5), "r", "LineWidth", 1.5)
 plot(x0 , y2 , "b--","LineWidth", 1.5)
 plot(x1 , y2 , "b-.","LineWidth", 1.5)
 plot(x2 , y2 , "r--","LineWidth", 1.5)
 plot(x3 , y2 , "r-.","LineWidth", 1.5)
 xlabel("Angulo de cigueñal[°]")
 ylabel("Densidad [kg/m^3]")
 title (" Densidad de mezcla Tubo 2 (adm)")
 print -dpdf densitube1d.pdf
 
 figure(2),%clf
 hold on, grid
 plot(tubo1a(:,3), tubo1a(:,6), "k", "LineWidth", 1.5)
 plot(tubo1b(:,3), tubo1b(:,6), "m", "LineWidth", 1.5)
 plot(tubo1c(:,3), tubo1c(:,6), "g", "LineWidth", 1.5)
 plot(x0 , y1 , "b--","LineWidth", 1.5)
 plot(x1 , y1 , "b-.","LineWidth", 1.5)
 plot(x2 , y1 , "r--","LineWidth", 1.5)
 plot(x3 , y1 , "r-.","LineWidth", 1.5)
 legend ("6000","7000", "8000","IVO","IVC","EVO","EVC") 
 xlabel("Angulo de cigueñal[°]")
 ylabel("Velocidad [m/seg]")
 title (" Velocidad de aire del tubo 2 ")
 print -dpdf prestube1d.pdf
 
 figure(3),%clf
 hold on, grid
 plot(tubo1a(:,3), tubo1a(:,7), "g", "LineWidth", 1.5)
 plot(x0 , y3 , "b--","LineWidth", 1.5)
 plot(x1 , y3 , "b-.","LineWidth", 1.5)
 plot(x2 , y3 , "r--","LineWidth", 1.5)
 plot(x3 , y3 , "r-.","LineWidth", 1.5)
 xlabel("Angulo de cigueñal[°]")
 ylabel("Presion [Pa]")
 title ("Presion del tubo 2")
 print -dpdf presitube1d.pdf
 
 