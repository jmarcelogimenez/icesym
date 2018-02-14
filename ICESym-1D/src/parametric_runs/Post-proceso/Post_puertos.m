%Script para el Post proceso de los puertos de adm y esc

clear
clc

rpms = [ 8000 ];

for irpm = 1:length(rpms)
  rpm = rpms(irpm);
  system(["grep '^  1  6' RPM_" num2str(rpm) "/cyl_0.txt > admision.txt"]);    
  adm = load(["admision.txt"]);
  
  system(["grep '^  2  6' RPM_" num2str(rpm) "/cyl_0.txt > escape.txt"]);    
  esc = load(["escape.txt"]);
  end
   
 figure(1),%clf
 hold on, grid
 plot(adm(:,3), adm(:,6), "r", "LineWidth", 1.5)
 plot(esc(:,3), esc(:,6), "b", "LineWidth", 1.5)
 xlabel("Angulo de cigueñal[°]")
 ylabel("Velocidad [m/s]")
 legend({"Admision","Escape"})
 title("Velocidad en puertos")
 print -dpdf velopuer.pdf
 
  figure(2),%clf
 hold on, grid
 plot(adm(:,3), adm(:,7), "r", "LineWidth", 1.5)
 plot(esc(:,3), esc(:,7), "b", "LineWidth", 1.5)
 xlabel("Angulo de cigueñal[°]")
 ylabel("Presion [Pa]")
 legend({"Admision","Escape"})
 title("Presion en puertos")
 print -dpdf uer.pdf