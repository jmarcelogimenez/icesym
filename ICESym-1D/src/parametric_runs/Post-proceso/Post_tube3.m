%Script para el Post proceso del tubo 3 en X=0

clear
clc

rpms = [ 6000 7000 8000];

  rpm = rpms(1);
  system(["grep '^  0  3' RPM_" num2str(rpm) "/tube_3.txt > tuboesc1.txt"]);    
  tubo3a = load(["tuboesc1.txt"]);
  
  rpm = rpms(2);
  system(["grep '^  0  3' RPM_" num2str(rpm) "/tube_3.txt > tuboesc2.txt"]);    
  tubo3b = load(["tuboesc2.txt"]);
  
  rpm = rpms(3);
  system(["grep '^  0  3' RPM_" num2str(rpm) "/tube_3.txt > tuboesc3.txt"]);    
  tubo3c = load(["tuboesc3.txt"]);
   
  velocidad_max= max(tubo3c(:,6));
  velocidad_min= min(tubo3c(:,6));
  
  IVO = 11.81
  IVC = 4.4156
  EVO = 8.22
  EVC = 0.82
  
  x0=[IVO/pi*180,IVO/pi*180];
  x1=[IVC/pi*180,IVC/pi*180];
  x2=[EVO/pi*180,EVO/pi*180];
  x3=[EVC/pi*180,EVC/pi*180];
  y1=[velocidad_min,velocidad_max];
  
  x4=[0, 720]
  y4=[101300, 101300] 
 
 figure(1),%clf
 subplot (2,2,1)  
 hold on, grid
 plot(tubo3c(:,3), tubo3c(:,5), "r", "LineWidth", 1.5)
 xlabel("Angulo de cigueñal[°]")
 ylabel("Densidad [kg/m^3]")
 title ("Densidad del tubo 3")
 
 subplot (2,2,2)
 hold on, grid
 plot(tubo3a(:,3), tubo3a(:,6), "k", "LineWidth", 1.5)
 plot(tubo3b(:,3), tubo3b(:,6), "m", "LineWidth", 1.5)
 plot(tubo3c(:,3), tubo3c(:,6), "g", "LineWidth", 1.5)
 plot(x0 , y1 , "b--","LineWidth", 1.5)
 plot(x1 , y1 , "b-.","LineWidth", 1.5)
 plot(x2 , y1 , "r--","LineWidth", 1.5)
 plot(x3 , y1 , "r-.","LineWidth", 1.5)
 legend ("6000","7000", "8000")
 xlabel("Angulo de cigueñal[°]")
 ylabel("Velocidad [m/seg]")
 title ("Velocidad del aire del tubo 3 en X=0")
 
subplot (2,2,3)
 hold on, grid
 plot(tubo3a(:,3), tubo3a(:,7), "k", "LineWidth", 1.5)
 plot(tubo3b(:,3), tubo3b(:,7), "m", "LineWidth", 1.5)
 plot(tubo3c(:,3), tubo3c(:,7), "g", "LineWidth", 1.5)
 plot(x4 , y4 , "b--","LineWidth", 1.5)
 xlabel("Angulo de cigueñal[°]")
 ylabel("Presion [Pa]")
 title ("Presion del tubo 3 en X=0")
 print -dpdf tube3.pdf