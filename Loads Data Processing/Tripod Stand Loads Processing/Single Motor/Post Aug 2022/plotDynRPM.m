function [f1,f2,f3,f4] = plotDynRPM(DynRPM,loads,GR)
%{
EDITED ON: 09/20/2022
EDITED BY: MATT ASPER

DETAILS: This update plots estimated motor torque with IQ

%}
%This function plots dynamic rpm data from Post-Aug 2022 Tripod Stand Testing


%% Plot
load('colors.mat')
xlimits_time = [-1 3];
xlimits_rev = [-13 39];

%Speed
%versus time
f1 = figure('Name','Rotor Speed');
subplot(2,1,1)
hold on
plot_areaerrorbar(DynRPM.time,DynRPM.Speed_avg/GR,DynRPM.Speed_err/GR,colors{1})
hold off
ylabel('$\Omega$, RPM')
legend('Experiment','location','southeast')
grid on
grid minor
xlim(xlimits_time)
xlabel('Time, s')

%versus revolution
subplot(2,1,2)
hold on
plot_areaerrorbar(DynRPM.rev/360/GR,DynRPM.Speed_avg/GR,DynRPM.Speed_err/GR,colors{1})
ylabel('$\Omega$, RPM')
legend('Experiment','location','southeast')
hold off
grid on
grid minor
xlim(xlimits_rev)
xlabel('Rotor Revolution')

%Torque
%versus time
f2 = figure('Name','Rotor Torque');
subplot(2,1,1)
if loads == 'y'; plot_areaerrorbar(DynRPM.time,DynRPM.Torque_avg*GR,DynRPM.Torque_err*GR,colors{1}); hold on; end
plot_areaerrorbar(DynRPM.time,DynRPM.Q_est_avg*GR,DynRPM.Q_est_err*GR,colors{2})
ylabel('$Q_r$, N$\cdot$m')
xlabel('Time, s')
if loads == 'y'; legend('Experiment','Estimated','location','southeast'); else; legend('Estimated','location','southeast'); end
formatfig
grid on
grid minor
xlim(xlimits_time)

%versus revolution
subplot(2,1,2)
if loads == 'y'; plot_areaerrorbar(DynRPM.rev/360/GR,DynRPM.Torque_avg*GR,DynRPM.Torque_err*GR,colors{1}); hold on; end
plot_areaerrorbar(DynRPM.rev/360/GR,DynRPM.Q_est_avg*GR,DynRPM.Q_est_err*GR,colors{2})
hold off
ylabel('$Q_r$, N$\cdot$m')
if loads == 'y'; legend('Experiment','Estimated','location','southeast'); else; legend('Estimated','location','southeast'); end
formatfig
grid on
grid minor
xlim(xlimits_rev)
xlabel('Rotor Revolution')
f2.Position = [326,236,674,561];

%Thrust
%versus time
f3 = figure('Name','Thrust');
subplot(2,1,1)
plot_areaerrorbar(DynRPM.time,DynRPM.T_avg,DynRPM.T_err,colors{1})
ylabel('T, N')
xlabel('Time, s')
formatfig
grid on
grid minor
xlim(xlimits_time)

%versus revolution
subplot(2,1,2)
plot_areaerrorbar(DynRPM.rev/360/GR,DynRPM.T_avg,DynRPM.T_err,colors{1})
ylabel('T, N')
formatfig
grid on
grid minor
xlim(xlimits_rev)
xlabel('Rotor Revolution')
f3.Position = [326,236,674,561];


%Torques
%versus time
f4 = figure('Name','Torque');
subplot(3,1,1)
plot_areaerrorbar(DynRPM.time,DynRPM.Torque_avg*GR,DynRPM.Torque_err*GR,colors{1})
ylabel('$Q_r$, N$\cdot$m')
formatfig

subplot(3,1,2)
plot_areaerrorbar(DynRPM.time,DynRPM.Qaero_avg,DynRPM.Qaero_err,colors{1})
ylabel('$Q_{aero}$, N$\cdot$m')
formatfig

subplot(3,1,3)
plot_areaerrorbar(DynRPM.time,DynRPM.Torque_avg*GR - DynRPM.Qaero_avg,DynRPM.Torque_err*GR + DynRPM.Qaero_err,colors{1})
ylabel('$Q_{I\dot{\alpha}}$, N$\cdot$m')
formatfig

xlabel('Time, s')
grid on
grid minor
xlim(xlimits_time)
f4.Position = [326,236,674,561];

%Calc inertia
f5 = figure('Name','Inertia');
plot(DynRPM.time,DynRPM.Inertia,'-','color',colors{1})
ylabel('I, N$\cdot$m')
xlabel('Time, s')
formatfig
xlim(xlimits_time)
f4.Position = [326,236,674,561];


end

