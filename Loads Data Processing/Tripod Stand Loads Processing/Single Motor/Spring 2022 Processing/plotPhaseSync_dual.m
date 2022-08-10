function [f1,f2,f3,f4,f5] = plotPhaseSync_dual(PhaseSync,loads)
%{
EDITED ON: 01/12/2022
EDITED BY: MATT ASPER

DETAILS: This update plots estimated motor torque with IQ

%}
%This function plots phase sync data from 2021-2022 Tripod Stand Testing

%% Constants 
%inertia (measured) kg*m^2
I = 0.3578;

%% Plot
load('colors.mat')

%Subplots
f1 = figure('Name','Response_vs_time');
subplot(3,1,1)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.ref_ang_avg,PhaseSync.ref_ang_err,colors{1})
hold off
ylabel('$\Delta\psi$, deg')
% legend('Experiment','','Prediction','location','northeast')
grid on
grid minor
xlim([-.1 .5])

subplot(3,1,2)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.Speed_avg,PhaseSync.Speed_err,colors{1})
hold off
ylabel('$\Omega$, RPM')
% legend('Experiment','','Prediction','location','northeast')
grid on
grid minor
xlim([-.1 .5])

subplot(3,1,3)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.Torque_avg,PhaseSync.Torque_err,colors{1});
hold off
ylabel('Q, N$\cdot$m')
xlabel('Time, s')
% legend('Experiment','location','southeast')
formatfig
grid on
grid minor
xlim([-.1 .5])
f1.Position = [326,236,674,561];

%Subplots
f2 = figure('Name','Response_vs_revolution');
subplot(3,1,1)
hold on
plot_areaerrorbar(PhaseSync.rev/360,PhaseSync.ref_ang_avg,PhaseSync.ref_ang_err,colors{1})
ylabel('$\Delta\psi$, deg')
% legend('Experiment','','Prediction','location','northeast')
hold off
grid on
grid minor
xlim([-1 10])

subplot(3,1,2)
hold on
plot_areaerrorbar(PhaseSync.rev/360,PhaseSync.Speed_avg,PhaseSync.Speed_err,colors{1})
ylabel('$\Omega$, RPM')
% legend('Experiment','','Prediction','location','northeast')
hold off
grid on
grid minor
xlim([-1 10])

subplot(3,1,3)
hold on
plot_areaerrorbar(PhaseSync.rev/360,PhaseSync.Torque_avg,PhaseSync.Torque_err,colors{1})
hold off
ylabel('Q, N$\cdot$m')
xlabel('Rotor Revolution')
% legend('Experiment','location','southeast')
formatfig
grid on
grid minor
xlim([-1 10])
f2.Position = [326,236,674,561];


%Thrust vs time
f3 = figure('Name','Thrust versus Time');
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.T_avg,PhaseSync.T_err,colors{1})
hold off
xlabel('Time, s')
ylabel('Thrust, N')
formatfig
xlim([-.1 1])

%Cts vs time
f4 = figure('Name','Blade Loading versus Time');
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.Cts_avg,PhaseSync.Cts_err,colors{1})
hold off
xlabel('Time, s')
ylabel('$C_T/\sigma$')
formatfig
xlim([-.1 1])

%Inertia vs time
f5 = figure('Name','Inertial Torque versus Time');
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.accel_avg*I,PhaseSync.accel_err*I,colors{1})
hold off
xlabel('Time, s')
ylabel('Q, N$\cdot$m')
formatfig
xlim([-.1 1])

end

