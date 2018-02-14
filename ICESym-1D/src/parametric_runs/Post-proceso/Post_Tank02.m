%Script para el Post proceso del tanque2 (SILENCIADOR)
% Calculo de Temperaturas vs Angulo de cigueñal
% Calculo de temperaturas promedio vs RPM 

clear
clc

rpms = [ 6000 6200 6400 6600 6800 7000 7200 7400 7500 7600 7700 7800 7900 8000];

  for i=1:length(rpms)
  rpm = rpms(i);
  system(["grep '^  0  6' RPM_" num2str(rpm) "/tank_2.txt > tanque2.txt"]);    
  matriz= load(["tanque2.txt"]);
  datos(:,i)=matriz(:,7); %Matriz de (Angulos de cig/filas) y (rpm/columnas)
  promedios(i)=mean(datos(:,i));
  end
  
  punto = load(["puestapto.txt"]);
    
  IVO = punto(1);
  IVC = punto(2);
  EVO = punto(3);
  EVC = punto(4);
  
  x0=[IVO/pi*180,IVO/pi*180];
  x1=[IVC/pi*180,IVC/pi*180];
  x2=[EVO/pi*180,EVO/pi*180];
  x3=[EVC/pi*180,EVC/pi*180];
  y=[max(datos(:,i)),min(datos(:,1))];
  
  figure(1)
  subplot(2,1,1)
  hold on, grid
  plot(matriz(:,3), datos(:,1), "r", "LineWidth", 1.5)    %Temperatura a 6000 RPM
  plot(matriz(:,3), datos(:,6), "b", "LineWidth", 1.5)    %Temperatura a 7000 RPM
  plot(matriz(:,3), datos(:,14), "g", "LineWidth", 1.5)   %Temperatura a 8000 RPM
  plot(x0 , y , "b--","LineWidth", 1.5)
  plot(x1 , y , "b--","LineWidth", 1.5)
  plot(x2 , y , "r--","LineWidth", 1.5)
  plot(x3 , y , "r--","LineWidth", 1.5)
  legend("6000","7000","8000")
  xlabel("Angulo de cigueñal[°]")
  ylabel("Temperatura [°K]")
  title ("Temperatura en el silenciador")
  
  
  subplot(2,1,2)
  hold on, grid
  plot(rpms, promedios, "k", "LineWidth", 1)
  xlabel("RPM")
  ylabel("Temperatura [°K]")
  title ("Temperatura promedio en el silenciador")
  print -dpdf TempTank2.jpg