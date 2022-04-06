function [f1,f2,f3] = plotSteady(Averages,collective)
%{
EDITED ON; 03/16/2022
EDITED BY: MATT ASPER

DETAILS:  This update plots dual motor steady data
%}

%Run this function to plot steady tripod 2021-2022 data from Averages
%Struct after running runTripod.m code

load('colors.mat')

%% Load Prediction

%call BEMT (Ct,Cp,theta)
BEMT = load('BEMTdata_coax_final.mat');

%create predicted Thrust and Power curves
RPM = linspace(0,1200,2000);
Ct_BEMT = interp1(BEMT.theta,BEMT.Ct,collective);
Cp_BEMT = interp1(BEMT.theta,BEMT.Cp,collective);
T_BEMT = Ct_BEMT*Averages.rho{1}*pi*Averages.R^2*(RPM*2*pi/60*Averages.R).^2;
P_BEMT = Cp_BEMT*Averages.rho{1}*pi*Averages.R^2*(RPM*2*pi/60*Averages.R).^3;


%% Plotting
f1 = figure('Name','Steady Ct versus Index Angle');
hold on
yline(Ct_BEMT/Averages.sigma,'k--','Linewidth',1.5)
errorbar(cell2mat(Averages.index_avg),cell2mat(Averages.cts_up_avg),...
    cell2mat(Averages.cts_up_err),cell2mat(Averages.cts_up_err),...
    cell2mat(Averages.index_err),cell2mat(Averages.index_err),'^','color',colors{1},...
    'MarkerFaceColor',colors{1})
errorbar(cell2mat(Averages.index_avg),cell2mat(Averages.cts_lo_avg),...
    cell2mat(Averages.cts_lo_err),cell2mat(Averages.cts_lo_err),...
    cell2mat(Averages.index_err),cell2mat(Averages.index_err),'o','color',colors{2},...
    'MarkerFaceColor',colors{2})
errorbar(cell2mat(Averages.index_avg),cell2mat(Averages.cts_avg),...
    cell2mat(Averages.cts_err),cell2mat(Averages.cts_err),...
    cell2mat(Averages.index_err),cell2mat(Averages.index_err),'s','color',colors{3},...
    'MarkerFaceColor',colors{3})
xlabel('Index Angle, $\phi$ [deg]')
ylabel('Blade Loading, $C_T/\sigma$ [-]')
legend('BEMT','Upper','Lower','Total','location','northeast')
% ylim([-inf (Cp_BEMT/Averages.sigma+0.01)])
formatfig
hold off

f2 = figure('Name','Steady Cp versus Index Angle');
hold on
yline(Cp_BEMT/Averages.sigma,'k--','Linewidth',1.5)
errorbar(cell2mat(Averages.index_avg),cell2mat(Averages.cps_up_avg),...
    cell2mat(Averages.cps_up_err),cell2mat(Averages.cps_up_err),...
    cell2mat(Averages.index_err),cell2mat(Averages.index_err),'^','color',colors{1},...
    'MarkerFaceColor',colors{1})
errorbar(cell2mat(Averages.index_avg),cell2mat(Averages.cps_lo_avg),...
    cell2mat(Averages.cps_lo_err),cell2mat(Averages.cps_lo_err),...
    cell2mat(Averages.index_err),cell2mat(Averages.index_err),'o','color',colors{2},...
    'MarkerFaceColor',colors{2})
errorbar(cell2mat(Averages.index_avg),cell2mat(Averages.cps_avg),...
    cell2mat(Averages.cps_err),cell2mat(Averages.cps_err),...
    cell2mat(Averages.index_err),cell2mat(Averages.index_err),'s','color',colors{3},...
    'MarkerFaceColor',colors{3})
xlabel('Index Angle, $\phi$ [deg]')
ylabel('Power Loading, $C_P/\sigma$ [-]')
legend('BEMT','Upper','Lower','Total','location','northeast')
ylim([-inf (Cp_BEMT/Averages.sigma+0.001)])
formatfig
hold off

f3 = figure('Name','Estimated Cps versus Index Angle');
hold on
yline(Cp_BEMT/Averages.sigma,'k--','Linewidth',1.5)
errorbar(cell2mat(Averages.index_avg),cell2mat(Averages.cps_up_est_avg),...
    cell2mat(Averages.cps_up_est_err),cell2mat(Averages.cps_up_est_err),...
    cell2mat(Averages.index_err),cell2mat(Averages.index_err),'^','color',colors{1},...
    'MarkerFaceColor',colors{1})
errorbar(cell2mat(Averages.index_avg),cell2mat(Averages.cps_lo_est_avg),...
    cell2mat(Averages.cps_lo_est_err),cell2mat(Averages.cps_lo_est_err),...
    cell2mat(Averages.index_err),cell2mat(Averages.index_err),'o','color',colors{2},...
    'MarkerFaceColor',colors{2})
errorbar(cell2mat(Averages.index_avg),cell2mat(Averages.cps_est_avg),...
    cell2mat(Averages.cps_est_err),cell2mat(Averages.cps_est_err),...
    cell2mat(Averages.index_err),cell2mat(Averages.index_err),'s','color',colors{3},...
    'MarkerFaceColor',colors{3})
xlabel('Index Angle, $\phi$ [deg]')
ylabel('Power Loading, $C_P/\sigma$ [-]')
legend('BEMT','Upper','Lower','Total','location','northeast')
ylim([-inf (Cp_BEMT/Averages.sigma+0.003)])
formatfig
hold off


end

