clc
close all
gpath = "";

B = 0.0855;       % bore diameter
L = 0.088;        % stroke carrera
rc = 9.5;       % compression ratio
rho_atm = 1.1647; % atmospheric density

Q_HV = 46.0e6;

n_R = 4;

Vd = pi/4*B^2*L;

rpms = [6000 6200 6400 6600 6800 7000 7200 7400 7500 7600 7700 7800 7900 8000 8100 8200];
lc   = 6; % last cycle
nac  = 1; % number of averaged cycles, from the last cycle   

IVC = 4.41568*180/pi;



if (exist('rpm', "var") == 0 || rpm<rpms(end))

for irpm = 1:length(rpms)

  rpm = rpms(irpm);

  Sp = 2*L*rpm/60;

  filename = [gpath "RPM_" num2str(rpm) "/cyl_extras_0.txt"];
  [theta0, time0, hc0, hr0, dQhc0, dQhr0, xb0, xbdot0, m_in0, m_out0, ...
   vol0, vdot0, m_fuel0, m_air0, m_res0, dQ_ht0, dQ_chem0] = ...
  read_cylextra(filename);
%  filename = [gpath "RPM_" num2str(rpm) "/cyl_extras_1.txt"];
%  [theta1, time1, hc1, hr1, dQhc1, dQhr1, xb1, xbdot1, m_in1, m_out1, ...
%   vol1, vdot1, m_fuel1, m_air1, m_res1, dQ_ht1, dQ_chem1] = ...
%  read_cylextra(filename);
%  filename = [gpath "RPM_" num2str(rpm) "/cyl_extras_2.txt"];
%  [theta2, time2, hc2, hr2, dQhc2, dQhr2, xb2, xbdot2, m_in2, m_out2, ...
%   vol2, vdot2, m_fuel2, m_air2, m_res2, dQ_ht2, dQ_chem2] = ...
%  read_cylextra(filename);
%  filename = [gpath "RPM_" num2str(rpm) "/cyl_extras_3.txt"];
%  [theta3, time3, hc3, hr3, dQhc3, dQhr3, xb3, xbdot3, m_in3, m_out3, ...
%   vol3, vdot3, m_fuel3, m_air3, m_res3, dQ_ht3, dQ_chem3] = ...
%  read_cylextra(filename);

  cyl0_state = load([gpath "RPM_" num2str(rpm) "/cyl_0.txt"]);
%  cyl1_state = load([gpath "RPM_" num2str(rpm) "/cyl_1.txt"]);
%  cyl2_state = load([gpath "RPM_" num2str(rpm) "/cyl_2.txt"]);
%  cyl3_state = load([gpath "RPM_" num2str(rpm) "/cyl_3.txt"]);

  cycle_work(irpm,4) = 0;
  mf_tr(irpm,4) = 0;
  ma_tr(irpm,4) = 0;
  mr_tr(irpm,4) = 0;
  Qht(irpm,4)   = 0;

  pcyl_max(irpm,4) = 0;

  rho_port(irpm,4) = 0;

  for icycle = lc-nac+1:lc
    ilc0 = find(cyl0_state(:,1)==0 & cyl0_state(:,2)==icycle);
%    ilc1 = find(cyl1_state(:,1)==0 & cyl1_state(:,2)==icycle);
%    ilc2 = find(cyl2_state(:,1)==0 & cyl2_state(:,2)==icycle);
%    ilc3 = find(cyl3_state(:,1)==0 & cyl3_state(:,2)==icycle);
    
    cycle_work(irpm,1) += trapz(vol0(1:length(ilc0),icycle),cyl0_state(ilc0,6));
%    cycle_work(irpm,2) += trapz(vol1(1:length(ilc1),icycle),cyl1_state(ilc1,6));
%    cycle_work(irpm,3) += trapz(vol2(1:length(ilc2),icycle),cyl2_state(ilc2,6));
%    cycle_work(irpm,4) += trapz(vol3(1:length(ilc3),icycle),cyl3_state(ilc3,6));

    pcyl_max(irpm,1) += max(cyl0_state(ilc0,6));
%    pcyl_max(irpm,2) += max(cyl1_state(ilc1,6));
%    pcyl_max(irpm,3) += max(cyl2_state(ilc2,6));
%    pcyl_max(irpm,4) += max(cyl3_state(ilc3,6));
    
    rho_port(irpm,1) += trapz(cyl0_state(ilc0,4),cyl0_state(ilc0,5))/(cyl0_state(ilc0(end),4)-cyl0_state(ilc0(1),4));
%    rho_port(irpm,2) += trapz(cyl1_state(ilc1,4),cyl1_state(ilc1,5))/(cyl1_state(ilc1(end),4)-cyl1_state(ilc1(1),4));
%    rho_port(irpm,3) += trapz(cyl2_state(ilc2,4),cyl2_state(ilc2,5))/(cyl2_state(ilc2(end),4)-cyl2_state(ilc2(1),4));
%    rho_port(irpm,4) += trapz(cyl3_state(ilc3,4),cyl2_state(ilc3,5))/(cyl2_state(ilc3(end),4)-cyl2_state(ilc3(1),4));  
      
    [dtmin0, imin0] = min(abs(theta0(:,icycle)-(IVC+1)));
    if(theta0(imin0,icycle) < IVC+1)
      imin0 += 1;
    end
    mf_tr(irpm,1) += m_fuel0(imin0,icycle);
    ma_tr(irpm,1) += m_air0(imin0,icycle);
    mr_tr(irpm,1) += m_res0(imin0,icycle);

%    [dtmin1, imin1] = min(abs(theta1(:,icycle)-(IVC+1)));
%    if(theta1(imin1,icycle) < IVC+1)
%      imin1 += 1;
%    end
%    mf_tr(irpm,2) += m_fuel1(imin1,icycle);
%    ma_tr(irpm,2) += m_air1(imin1,icycle);
%    mr_tr(irpm,2) += m_res1(imin1,icycle);
%
%    	[dtmin2, imin2] = min(abs(theta2(:,icycle)-(IVC+1)));
%    if(theta2(imin2,icycle) < IVC+1)
%      imin2 += 1;
%    end
%    mf_tr(irpm,3) += m_fuel2(imin2,icycle);
%    ma_tr(irpm,3) += m_air2(imin2,icycle);
%    mr_tr(irpm,3) += m_res2(imin2,icycle);
%
%    [dtmin3, imin3] = min(abs(theta3(:,icycle)-(IVC+1)));
%    if(theta3(imin3,icycle) < IVC+1)
%      imin3 += 1;
%    end
%    mf_tr(irpm,4) += m_fuel3(imin3,icycle);
%    ma_tr(irpm,4) += m_air3(imin3,icycle);
%    mr_tr(irpm,4) += m_res3(imin3,icycle);
    
    Qht(irpm,1) += trapz(theta0(:,icycle),dQ_ht0(:,icycle))/720;
%    Qht(irpm,2) += trapz(theta1(:,icycle),dQ_ht1(:,icycle))/720;
%    Qht(irpm,3) += trapz(theta2(:,icycle),dQ_ht2(:,icycle))/720;
%    Qht(irpm,4) += trapz(theta3(:,icycle),dQ_ht3(:,icycle))/720;
  end
  cycle_work(irpm,:) /= nac;
  pcyl_max(irpm,:) /= nac;
  # mf_tr(irpm,:) /= nac;
  mf_tr(irpm,:) = 5.1584e-5;
  ma_tr(irpm,:) /= nac;
  mr_tr(irpm,:) /= nac;
  Qht(irpm,:) /= nac;

  ipower(irpm,:)  = cycle_work(irpm,:)*rpm/60/n_R;
  itorque(irpm,:) = cycle_work(irpm,:)/(2*pi*n_R);
  imep(irpm,:)    = cycle_work(irpm,:)/Vd;

  ifuelconv_eff(irpm,:) = cycle_work(irpm,:)./(mf_tr(irpm,:)*Q_HV);
  isfc(irpm,:)         = 3600*1e6./(ifuelconv_eff(irpm,:)*Q_HV);

  fmep(irpm,:) = 1e5*(0.4+0.005e-5*pcyl_max(irpm,:)+0.09*Sp+9e-4*Sp^2);

  bmep(irpm,:)    = imep(irpm,:)-fmep(irpm,:);
  bpower(irpm,:)  = bmep(irpm,:)*Vd*rpm/60/n_R;
  btorque(irpm,:) = bmep(irpm,:)*Vd/(2*pi*n_R);

  bfuelconv_eff(irpm,:) = bmep(irpm,:)*Vd./(mf_tr(irpm,:)*Q_HV);
  bsfc(irpm,:)         = 3600*1e6./(bfuelconv_eff(irpm,:)*Q_HV);

  vol_eff(irpm,:)  = ma_tr(irpm,:)/(rho_atm*Vd);
  vol_eff1(irpm,:) = ma_tr(irpm,:)./(rho_port(irpm,:)*Vd);
  xr(irpm,:)       = mr_tr(irpm,:)./(mf_tr(irpm,:)+ma_tr(irpm,:)+mr_tr(irpm,:));

end

endif

ipower_total = sum(ipower)*0.00134102209; # HP
bpower_total = sum(bpower)*0.00134102209; # HP
isfc_total   = sum(isfc.*ipower)/sum(ipower);
bsfc_total   = sum(bsfc.*bpower)/sum(bpower);



figure(1),%clf
 hold on, grid
 plot(rpms,itorque(:,end),"b")
 plot(rpms,ipower(:,end)/735.5, "r")
 xlabel("engine speed [rpm]")
 ylabel("power [CV], torque [Nm]")
 legend("ipower", "itorque")
%
 figure(2),clf
 plot(rpms,vol_eff(:,end-1))
 hold on, grid
 plot(rpms,ifuelconv_eff(:,end-1),'r')
 plot(rpms,bfuelconv_eff(:,end-1),'g')
 plot(rpms,xr(:,end-1),"k")
 xlabel("engine speed [rpm]")
 legend("volumetric eff.", "i fuel conversion eff.","b fuel conversion eff.","residual mass fraction")

%figure(3),clf
%plot(rpms,isfc)
%grid
%xlabel("engine speed [rpm]")
%ylabel("specific fuel consumption [g/kW h]")
%legend("cycle 1", "cycle 2", "cycle 3", "cycle 4")
%
%figure(4),clf
%plot(rpms,imep/1000)
%grid
%xlabel("engine speed [rpm]")
%ylabel("mean effective pressure [kPa]")

figure(3),clf
plot(theta0(:,lc),m_in0(:,lc),"b")
hold on
plot(theta1(:,lc),m_in1(:,lc),"g")
plot(theta0(:,lc),-m_out0(:,lc),"r")
plot(theta1(:,lc),-m_out1(:,lc),"m")
grid
axis([0 720 -0.02 0.1])
xlabel("angulo del ciguenial [grad]")
ylabel("caudal masico [kg/s]")
%
figure(4),clf
plot(theta0(:,lc),m_in0(:,lc),"b")
hold on
plot(theta1(:,lc),m_in1(:,lc),"g")
plot(theta0(:,lc),-m_out0(:,lc),"r")
plot(theta1(:,lc),-m_out1(:,lc),"m")
grid
axis([0 720 -0.02 0.1])
xlabel("angulo del ciguenial [grad]")
ylabel("caudal masico [kg/s]")

%figure(2),clf
%plot(theta0(:,lc),xbdot0(:,lc),"b")
%hold on
%plot([339 339],[0 600])
%grid
%axis([300 500 0 600])
%xlabel("angulo del ciguenial [grad]")
%ylabel("dxb / dt [1/s]")
