function [f1,f2] = plotDynRPM(DynRPM,loads,GR)
%{
EDITED ON: 09/20/2022
EDITED BY: MATT ASPER

DETAILS: This update plots estimated motor torque with IQ

%}
%This function plots dynamic rpm data from Post-Aug 2022 Tripod Stand Testing


%% Plot
load('colors.mat')
xlimits_time = [-1 2];
xlimits_rev = [-13 26];

%Speed
%versus time
f1 = figure('Name','Rotor Speed');
subplot(2,1,1)
hold on
plot_areaerrorbar(DynRPM.time,DynRPM.Speed_avg/GR,DynRPM.Speed_err/GR,colors{1})
hold off
ylabel('$\Omega$, RPM')
legend('Experiment','','Prediction','location','southeast')
grid on
grid minor
xlim(xlimits_time)
xlabel('Time, s')

%versus revolution
subplot(2,1,2)
hold on
plot_areaerrorbar(DynRPM.rev/360/GR,DynRPM.Speed_avg/GR,DynRPM.Speed_err/GR,colors{1})
ylabel('$\Omega$, RPM')
legend('Experiment','','Prediction','location','southeast')
hold off
grid on
grid minor
xlim(xlimits_rev)
xlabel('Rotor Revolution')

%Torque
%versus time
f2 = figure('Name','Rotor Torque');
subplot(2,1,1)
hold on
if loads == 'y'; plot_areaerrorbar(DynRPM.time,DynRPM.Torque_avg*GR,DynRPM.Torque_err*GR,colors{1}); end
plot_areaerrorbar(DynRPM.time,DynRPM.Q_est_avg*GR,DynRPM.Q_est_err*GR,colors{2})
hold off
ylabel('$Q_r$, $N\cdot m$')
xlabel('Time, s')
if loads == 'y'; legend('Experiment','Estimated','location','southeast'); else; legend('Estimated','location','southeast'); end
formatfig
grid on
grid minor
xlim(xlimits_time)
f1.Position = [326,236,674,561];

%versus revolution
subplot(2,1,2)
hold on
if loads == 'y'; plot_areaerrorbar(DynRPM.rev/360/GR,DynRPM.Torque_avg*GR,DynRPM.Torque_err*GR,colors{1}); end
plot_areaerrorbar(DynRPM.rev/360/GR,DynRPM.Q_est_avg*GR,DynRPM.Q_est_err*GR,colors{2})
hold off
ylabel('$Q_r$, $N\cdot m$')
if loads == 'y'; legend('Experiment','Estimated','location','southeast'); else; legend('Estimated','location','southeast'); end
formatfig
grid on
grid minor
xlim(xlimits_rev)
xlabel('Rotor Revolution')
f2.Position = [326,236,674,561];



end

