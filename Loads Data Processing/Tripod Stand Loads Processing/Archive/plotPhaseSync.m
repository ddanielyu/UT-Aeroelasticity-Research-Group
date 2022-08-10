function [f1,f2,f3] = plotPhaseSync(PhaseSync,PhaseSync_test)
%{
EDITED ON: 01/12/2022
EDITED BY: MATT ASPER

DETAILS: This update plots estimated motor torque with IQ

%}
%This function plots phase sync data from 2021-2022 Tripod Stand Testing

%% Extract unique test ID's
for i = 1:length(PhaseSync_test)
    test_name = split(PhaseSync_test{i},{'_','-'});
    tests{i} = test_name{3};
end
unique_tests = unique(tests);

fprintf('\nPhase-sync tests: ')
for i = 1:length(unique_tests)
    fprintf('%s ',unique_tests{i})
end
whattoplot = input('\nTest to plot: ','s');

id = find(cell2mat(unique_tests) == whattoplot);


%% Fit Current-Torque Line
ft = fittype('a*x');
[fitobj1,~,~,~] = fit(PhaseSync.Curr1_avg(1:end-1)'/sqrt(2),PhaseSync.Torque_pk_avg(1:end-1)',ft);
curr_xline = linspace(0,max(PhaseSync.Curr1_avg(1:end-1)/sqrt(2))+10,2000);
curr_yline = fitobj1.a*curr_xline;

%% Plot
load('colors.mat')

%Subplots
f1 = figure('Name','Response_vs_time');
subplot(3,1,1)
hold on
plot_areaerrorbar(PhaseSync.time{id},PhaseSync.Ang_err_avg{id},PhaseSync.Ang_err_err{id},colors{1})
% plot(Angle_err.time, Angle_err.data * -1,'k-','linewidth',1.5)
hold off
ylabel('Angle Error, deg')
legend('Experiment','','Prediction')
grid on
grid minor
xlim([-.1 1])

subplot(3,1,2)
hold on
plot_areaerrorbar(PhaseSync.time{id},PhaseSync.Speed_avg{id},PhaseSync.Speed_err{id},colors{1})
% plot(Motor_RPM.time, Motor_RPM.data,'k-','linewidth',1.5)
hold off
ylabel('$\Omega_{motor}$, RPM')
legend('Experiment','','Prediction')
grid on
grid minor
xlim([-.1 1])

subplot(3,1,3)
hold on
plot_areaerrorbar(PhaseSync.time{id},PhaseSync.Torque_avg{id},PhaseSync.Torque_err{id},colors{1})
plot_areaerrorbar(PhaseSync.time{id},PhaseSync.Q_est_avg{id},PhaseSync.Q_est_err{id},colors{2})
% plot(Q_total.time, Q_total.data/GR,'k-','linewidth',1.5)
hold off
ylabel('$Q_{motor}$, $N\cdot m$')
xlabel('Time, s')
legend('Experiment','Estimated')
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
plot_areaerrorbar(PhaseSync.rev{id}/360,PhaseSync.Ang_err_avg{id},PhaseSync.Ang_err_err{id},colors{1})
% plot(nRev_new(full_step1:end)-nRev_new(full_step),Angle_err.data(full_step1:end) * -1,'k-','linewidth',1.5)
ylabel('Angle Error, deg')
legend('Experiment','','Prediction')
hold off
grid on
grid minor
xlim([-1 8])

subplot(3,1,2)
hold on
plot_areaerrorbar(PhaseSync.rev{id}/360,PhaseSync.Speed_avg{id},PhaseSync.Speed_err{id},colors{1})
% plot(nRev_new(full_step1:end)-nRev_new(full_step),Motor_RPM.data(full_step1:end-6),'k-','linewidth',1.5)
ylabel('$\Omega_{motor}$, RPM')
legend('Experiment','','Prediction')
hold off
grid on
grid minor
xlim([-1 8])

subplot(3,1,3)
hold on
plot_areaerrorbar(PhaseSync.rev{id}/360,PhaseSync.Torque_avg{id},PhaseSync.Torque_err{id},colors{1})
plot_areaerrorbar(PhaseSync.rev{id}/360,PhaseSync.Q_est_avg{id},PhaseSync.Q_est_err{id},colors{2})
% plot(nRev_new(full_step1:end)-nRev_new(full_step),Q_total.data(full_step1:end-6)/GR,'k-','linewidth',1.5)
hold off
ylabel('$Q_{motor}$, $N\cdot m$')
xlabel('Motor Revolution')
sgtitle('$+5^\circ$ Angle Offset','Fontsize',24)
legend('Experiment','Estimated')
formatfig
grid on
grid minor
xlim([-1 8])
f2.Position = [326,236,674,561];

f3 = figure('Name','Peak Currents versus Torque');
hold on
errorbar(PhaseSync.Curr1_avg(1:end-1)/sqrt(2),PhaseSync.Torque_pk_avg(1:end-1),PhaseSync.Torque_pk_err(1:end-1),PhaseSync.Torque_pk_err(1:end-1),...
    PhaseSync.Curr1_err(1:end-1)/sqrt(2),PhaseSync.Curr1_err(1:end-1)/sqrt(2),'^','color',colors{1})
plot(curr_xline,curr_yline,'k--')
curr_ax = gca;
text(min(PhaseSync.Curr1_avg(1:end-1)/sqrt(2)),max(PhaseSync.Torque_pk_avg(1:end-1)),strcat('Y = ',num2str(fitobj1.a),'X'),...
    'fontsize',16);
xlabel('Current, $A_{rms}$')
ylabel('Motor Torque, N$\cdot$m')
xlim([min(PhaseSync.Curr1_avg(1:end-1)/sqrt(2))-1,max(PhaseSync.Curr1_avg(1:end-1)/sqrt(2))+1])
formatfig
grid on
grid minor
hold off


end

