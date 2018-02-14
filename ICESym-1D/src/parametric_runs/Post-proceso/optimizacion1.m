%Scrip de optimizacion de largo de escape primario
clear clc

  B = 0.0865; #m
  a = 0.044; #m
  l = 0.14;	#m
  Vc = 0.000058; #m^3
  pi = 3.14159;
  Vd=pi*B^2/4*2*a;
  
rpms = [6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7500, 7600,  7700, 7800, 7900, 8000];
powData = [215.4, 225.2, 239, 249.4, 255.5, 260.9, 266.9, 272.4, 275.8, 279.4, 282.4, 283.6, 282.6, 283.4];
torData = [25.7, 26, 26.7, 27.1, 26.9, 26.7, 26.6, 26.4, 26.3, 26.3, 26.3, 26, 25.6, 25.4];

theta0 = [ 5.6985, 5.69, 5.683, 5.67, 5.66, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65, 5.65]; %%Extraido mapa Boscarol

L=[0 0]
inc=0
n=5;
potencia(1)=0

input("ingreso el largo de escape primario en mm (350 es el mejor): ")
L(2)=ans/1000
input("ingreso el delta longitud para iterar en mm: ")
inc=ans/1000

for z=2:n
    L4 = L(z)	%Prolongacion de escape, original = 0.25
    system(["sed -i 's/L4 =.*/L4 = " num2str(L4) ";/g' BoscaMono.py"]);
    
    dtheta = 0.002*2*pi*rpms/60

    PHI = [1.25	1.25	1.252	1.255	1.26	1.29	1.335	1.385	1.395	1.405	1.41	1.42	1.435	1.455]; %1.48 1.52] %1.25+0.2/2200*[rpms-6000]	%Modelo lineal PHI(rpms)

    for irpm = 1:length(rpms)
	  system(["sed -i 's/rpms =.*/rpms = [" num2str(rpms(irpm)) "];/g' BoscaMono.py"]);
	  system(["sed -i 's/PHI =.*/PHI = " num2str(PHI(irpm)) ";/g' BoscaMono.py"]);
	  system(["sed -i 's/theta0 =.*/theta0 = " num2str(theta0(irpm)) ";/g' BoscaMono.py"]);
	  system(["sed -i 's/dtheta =.*/dtheta = " num2str(dtheta(irpm)) ";/g' BoscaMono.py"]);
	  system(["python monoRun.py >> Mono.log"]); 
    end

    for irpm = 1:length(rpms)

	  rpm = rpms(irpm); Sp = 4*a*rpm/60;
	  system(["grep '^  0  3' tests/BoscaMono/RPM_" num2str(rpm) "/cyl_0.txt > extractedCyl0.txt"]);    
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
    
    pot(z,irpm)=power(irpm);
    tor(z,irpm)=torque(irpm);
    end

    potencia(z)=sum(pot(z,:))/length(rpms)
   

      if potencia(z) > potencia(z-1)
      n=z+1
      %n=n+1
      L(n)=L(n-1)+inc
      end   
 end
 
clear
disp("La longitud ideal en mm del primario es:"), L(z-1)*1000
 
figure(1),%clf
hold on ,grid
plot(rpms,pot(z-1,:)/735.5, "k")    %Mejor resultado
plot(rpms,pot(z,:)/735.5, "k**")    %Ultimo resultado
plot(rpms,powData(1:length(rpms))/4, "r--")
plot(rpms,torData(1:length(rpms))*9.8/4, "b--")
xlabel("engine speed [rpm]")
ylabel("power [CV], torque [Nm]")
legend({"Largo optimo", "Ultima largo ensayado", "test power", "test torque"}, 'location', 'southeast')
print -dpdf powerr.pdf  
    