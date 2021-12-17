function plotPhaseSync(PhaseSync)
%This function plots phase sync data from 2021-2022 Tripod Stand Testing

load('colors.mat')

%Subplots
f1 = figure('Name','Response_vs_time');
subplot(3,1,1)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.Ang_err_avg,PhaseSync.Ang_err_err,colors{1})
% plot(Angle_err.time, Angle_err.data * -1,'k-','linewidth',1.5)
hold off
ylabel('Angle Error, deg')
legend('Experiment','','Prediction')
grid on
grid minor
xlim([-.1 1])

subplot(3,1,2)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.Speed_avg,PhaseSync.Speed_err,colors{1})
% plot(Motor_RPM.time, Motor_RPM.data,'k-','linewidth',1.5)
hold off
ylabel('$\Omega_{motor}$, RPM')
legend('Experiment','','Prediction')
grid on
grid minor
xlim([-.1 1])

subplot(3,1,3)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.Torque_avg,PhaseSync.Torque_err,colors{1})
% plot(Q_total.time, Q_total.data/GR,'k-','linewidth',1.5)
hold off
ylabel('$Q_{motor}$, $N\cdot m$')
xlabel('Time, s')
legend('Experiment','','Prediction')
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
plot_areaerrorbar(PhaseSync.rev/360,PhaseSync.Ang_err_avg,PhaseSync.Ang_err_err,colors{1})
% plot(nRev_new(full_step1:end)-nRev_new(full_step),Angle_err.data(full_step1:end) * -1,'k-','linewidth',1.5)
ylabel('Angle Error, deg')
legend('Experiment','','Prediction')
hold off
grid on
grid minor
xlim([-1 8])

subplot(3,1,2)
hold on
plot_areaerrorbar(PhaseSync.rev/360,PhaseSync.Speed_avg,PhaseSync.Speed_err,colors{1})
% plot(nRev_new(full_step1:end)-nRev_new(full_step),Motor_RPM.data(full_step1:end-6),'k-','linewidth',1.5)
ylabel('$\Omega_{motor}$, RPM')
legend('Experiment','','Prediction')
hold off
grid on
grid minor
xlim([-1 8])

subplot(3,1,3)
hold on
plot_areaerrorbar(PhaseSync.rev/360,PhaseSync.Torque_avg,PhaseSync.Torque_err,colors{1})
% plot(nRev_new(full_step1:end)-nRev_new(full_step),Q_total.data(full_step1:end-6)/GR,'k-','linewidth',1.5)
hold off
ylabel('$Q_{motor}$, $N\cdot m$')
xlabel('Motor Revolution')
sgtitle('$+5^\circ$ Angle Offset','Fontsize',24)
legend('Experiment','','Prediction')
formatfig
grid on
grid minor
xlim([-1 8])
f2.Position = [326,236,674,561];

end

