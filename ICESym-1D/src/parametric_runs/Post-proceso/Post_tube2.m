%Script para el Post proceso del tubo 2 (multiple de admision)

clear
clc

rpms = [ 6000 6200 6400 6600 6800 7000 7200 7400 7500 7600 7700 7800 7900 8000];

  for i=1:length(rpms)
  rpm = rpms(i);
  system(["grep '^  0  6' RPM_" num2str(rpm) "/tube_2.txt > tubo2.txt"]);    
  matriz = load(["tubo2.txt"]);
  datosPres(:,i)=matriz(:,7);
  datosVel(:,i)=matriz(:,6);
  datosDen(:,i)=matriz(:,5);
  end

  presion_max= max(datosPres(:,i));
  presion_min= min(datosPres(:,1));
  densidad_max= max(datosDen(:,i));
  densidad_min= min(datosDen(:,1));
  velocidad_max= max(datosVel(:,i));
  velocidad_min= min(datosVel(:,1));
  
  punto = load(["puestapto.txt"]);
    
  IVO = punto(1);
  IVC = punto(2);
  EVO = punto(3);
  EVC = punto(4);
  
  x0=[IVO/pi*180,IVO/pi*180];
  x1=[IVC/pi*180,IVC/pi*180];
  x2=[EVO/pi*180,EVO/pi*180];
  x3=[EVC/pi*180,EVC/pi*180];
  y1=[velocidad_min,velocidad_max];
  y2=[densidad_min,densidad_max];
  y3=[presion_min,presion_max];  
   
 figure(1),%clf
 hold on, grid
 plot(matriz(:,3), datosDen(:,1), "r", "LineWidth", 1)  %Densidad a 6000 RPM
 plot(matriz(:,3), datosDen(:,6), "b", "LineWidth", 1)  %Densidad a 7000 RPM
 plot(matriz(:,3), datosDen(:,14), "g", "LineWidth", 1) %Densidad a 8000 RPM
 plot(x0 , y2 , "b--","LineWidth", 1.5)
 plot(x1 , y2 , "b--","LineWidth", 1.5)
 plot(x2 , y2 , "r--","LineWidth", 1.5)
 plot(x3 , y2 , "r--","LineWidth", 1.5)
 legend("6000","7000","8000");
 xlabel("Angulo de cigueñal[°]")
 ylabel("Densidad [kg/m^3]")
 title (" Densidad de mezcla Tubo 2 (adm)")
 print -dpdf Tube2-Dens.pdf
 
 figure(2),%clf
 hold on, grid
 plot(matriz(:,3), datosVel(:,1), "r", "LineWidth", 1)  %Velocidad a 6000 RPM
 plot(matriz(:,3), datosVel(:,6), "b", "LineWidth", 1)  %Velocidad a 7000 RPM
 plot(matriz(:,3), datosVel(:,14), "g", "LineWidth", 1) %Velocidad a 8000 RPM
 plot(x0 , y1 , "b--","LineWidth", 1.5)
 plot(x1 , y1 , "b--","LineWidth", 1.5)
 plot(x2 , y1 , "r--","LineWidth", 1.5)
 plot(x3 , y1 , "r--","LineWidth", 1.5)
 legend("6000","7000","8000");
 xlabel("Angulo de cigueñal[°]")
 ylabel("Velocidad [m/seg]")
 title (" Velocidad de aire en Tubo 2 (adm)")
 print -dpdf Tube2-Vel.pdf

 
 figure(3),%clf
 hold on, grid
 plot(matriz(:,3), datosPres(:,1), "r", "LineWidth", 1)  %Presion a 6000 RPM
 plot(matriz(:,3), datosPres(:,6), "b", "LineWidth", 1)  %Presion a 7000 RPM
 plot(matriz(:,3), datosPres(:,14), "g", "LineWidth", 1) %Presion a 8000 RPM
 plot(x0 , y3 , "b--","LineWidth", 1.5)
 plot(x1 , y3 , "b--","LineWidth", 1.5)
 plot(x2 , y3 , "r--","LineWidth", 1.5)
 plot(x3 , y3 , "r--","LineWidth", 1.5)
 legend("6000","7000","8000");
 xlabel("Angulo de cigueñal[°]")
 ylabel("Presion [Pa]")
 title (" Presion del aire en Tubo 2 (adm)")
 print -dpdf Tube2-Pres.pdf %depsc para guardar como imagen en color
 
