function [f1,f2,f3] = plotPhaseSync_dual(PhaseSync,loads)
%{
EDITED ON: 01/12/2022
EDITED BY: MATT ASPER

DETAILS: This update plots estimated motor torque with IQ

%}
%This function plots phase sync data from 2021-2022 Tripod Stand Testing

%% Plot
load('colors.mat')

%Subplots
f1 = figure('Name','Response_vs_time');
subplot(3,1,1)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.ref_ang_avg,PhaseSync.ref_ang_err,colors{1})
hold off
ylabel('Motor Angle, deg')
legend('Experiment','','Prediction','location','southeast')
grid on
grid minor
xlim([-.1 1])

subplot(3,1,2)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.Speed_avg,PhaseSync.Speed_err,colors{1})
hold off
ylabel('$\Omega_{rotor}$, RPM')
legend('Experiment','','Prediction','location','northeast')
grid on
grid minor
xlim([-.1 1])

subplot(3,1,3)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.Torque_avg,PhaseSync.Torque_err,colors{1});
hold off
ylabel('$Q_{rotor}$, $N\cdot m$')
xlabel('Time, s')
legend('Experiment','location','southeast')
sgtitle('$+5^\circ$ Angle Offset','Fontsize',24)
formatfig
grid on
grid minor
xlim([-.1 1])
f1.Position = [326,236,674,561];

%Subplots
f2 = figure('Name','Response_vs_revolution');
subplot(3,1,1)
hold on
plot_areaerrorbar(PhaseSync.rev/360,PhaseSync.ref_ang_avg,PhaseSync.ref_ang_err,colors{1})
ylabel('Motor Angle, deg')
legend('Experiment','','Prediction','location','southeast')
hold off
grid on
grid minor
xlim([-1 10])

subplot(3,1,2)
hold on
plot_areaerrorbar(PhaseSync.rev/360,PhaseSync.Speed_avg,PhaseSync.Speed_err,colors{1})
ylabel('$\Omega_{rotor}$, RPM')
legend('Experiment','','Prediction','location','northeast')
hold off
grid on
grid minor
xlim([-1 10])

subplot(3,1,3)
hold on
plot_areaerrorbar(PhaseSync.rev/360,PhaseSync.Torque_avg,PhaseSync.Torque_err,colors{1})
hold off
ylabel('$Q_{rotor}$, $N\cdot m$')
xlabel('Motor Revolution')
legend('Experiment','location','southeast')
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

end

