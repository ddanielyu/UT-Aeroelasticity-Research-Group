function plotSteady(Averages,collective)
%Run this function to plot steady tripod 2021-2022 data from Averages
%Struct after running runTripod.m code

load('colors.mat')

%% Load Prediction

%call BEMT (Ct,Cp,theta)
BEMT = load('/Users/asper101/Box Sync/Matt Lab Stuff/3rd Year/BEMT/FINAL/BEMTdata_rev2.mat');

%create predicted Thrust and Power curves
RPM = linspace(0,1200,2000);
Ct_BEMT = interp1(BEMT.theta,BEMT.Ct,collective);
Cp_BEMT = interp1(BEMT.theta,BEMT.Cp,collective);
T_BEMT = Ct_BEMT*Averages.rho{1}*pi*Averages.R^2*(RPM*2*pi/60*Averages.R).^2;
P_BEMT = Cp_BEMT*Averages.rho{1}*pi*Averages.R^2*(RPM*2*pi/60*Averages.R).^3;


%% Plotting
figure('Name','Steady Ct')
hold on
plot(BEMT.theta,BEMT.Ct,'k-','Linewidth',1.5)
errorbar(ones(1,length(Averages.cts_avg))*collective,cell2mat(Averages.cts_avg)*Averages.sigma,cell2mat(Averages.cts_err)*Averages.sigma,'^','color',colors{1})
xlabel('Collective, $\theta_o$ [deg]')
ylabel('Thrust Coefficient, $C_T$ [-]')
formatfig
hold off

figure('Name','Steady Cp')
hold on
plot(BEMT.theta,BEMT.Cp,'k-','Linewidth',1.5)
errorbar(ones(1,length(Averages.cps_avg))*collective,cell2mat(Averages.cps_avg)*Averages.sigma,cell2mat(Averages.cps_err)*Averages.sigma,'^','color',colors{1})
xlabel('Collective, $\theta_o$ [deg]')
ylabel('Power Coefficient, $C_P$ [-]')
formatfig
hold off

figure('Name','Steady Thrust')
hold on
plot(RPM,T_BEMT,'k-','Linewidth',1.5)
errorbar(cell2mat(Averages.OMEGA)*60/2/pi,cell2mat(Averages.T_avg),cell2mat(Averages.T_err),'^','color',colors{1})
xlabel('Rotor Speed, $\Omega$ [RPM]')
ylabel('Thrust [N]')
xlim([min(cell2mat(Averages.OMEGA)*60/2/pi) max(cell2mat(Averages.OMEGA)*60/2/pi)])
formatfig
hold off

figure('Name','Steady Power')
hold on
plot(RPM,P_BEMT,'k-','Linewidth',1.5)
errorbar(cell2mat(Averages.OMEGA)*60/2/pi,cell2mat(Averages.P_avg),cell2mat(Averages.P_err),'^','color',colors{1})
xlabel('Rotor Speed, $\Omega$ [RPM]')
ylabel('Power [W]')
xlim([min(cell2mat(Averages.OMEGA)*60/2/pi) max(cell2mat(Averages.OMEGA)*60/2/pi)])
formatfig
hold off

end

