function [f1,f2] = plotPhaseSync(PhaseSync,PhaseSync_test,loads)
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

%% Plot
load('colors.mat')

%Subplots
f1 = figure('Name','Response_vs_time');
subplot(3,1,1)
hold on
plot_areaerrorbar(PhaseSync.time{id},PhaseSync.ref_ang_avg{id},PhaseSync.ref_ang_err{id},colors{1})
hold off
ylabel('Motor Angle, deg')
legend('Experiment','','Prediction','location','southeast')
grid on
grid minor
xlim([-.1 1])

subplot(3,1,2)
hold on
plot_areaerrorbar(PhaseSync.time{id},PhaseSync.Speed_avg{id},PhaseSync.Speed_err{id},colors{1})
hold off
ylabel('$\Omega_{motor}$, RPM')
legend('Experiment','','Prediction','location','northeast')
grid on
grid minor
xlim([-.1 1])

subplot(3,1,3)
hold on
if loads == 'y'; plot_areaerrorbar(PhaseSync.time{id},PhaseSync.Torque_avg{id},PhaseSync.Torque_err{id},colors{1}); end
plot_areaerrorbar(PhaseSync.time{id},PhaseSync.Q_est_avg{id},PhaseSync.Q_est_err{id},colors{2})
hold off
ylabel('$Q_{motor}$, $N\cdot m$')
xlabel('Time, s')
if loads == 'y'; legend('Experiment','Estimated','location','southeast'); else; legend('Estimated','location','southeast'); end
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
plot_areaerrorbar(PhaseSync.rev{id}/360,PhaseSync.ref_ang_avg{id},PhaseSync.ref_ang_err{id},colors{1})
ylabel('Motor Angle, deg')
legend('Experiment','','Prediction','location','southeast')
hold off
grid on
grid minor
xlim([-1 10])

subplot(3,1,2)
hold on
plot_areaerrorbar(PhaseSync.rev{id}/360,PhaseSync.Speed_avg{id},PhaseSync.Speed_err{id},colors{1})
ylabel('$\Omega_{motor}$, RPM')
legend('Experiment','','Prediction','location','northeast')
hold off
grid on
grid minor
xlim([-1 10])

subplot(3,1,3)
hold on
if loads == 'y'; plot_areaerrorbar(PhaseSync.rev{id}/360,PhaseSync.Torque_avg{id},PhaseSync.Torque_err{id},colors{1}); end
plot_areaerrorbar(PhaseSync.rev{id}/360,PhaseSync.Q_est_avg{id},PhaseSync.Q_est_err{id},colors{2})
hold off
ylabel('$Q_{motor}$, $N\cdot m$')
xlabel('Motor Revolution')
sgtitle('$+5^\circ$ Angle Offset','Fontsize',24)
if loads == 'y'; legend('Experiment','Estimated','location','southeast'); else; legend('Estimated','location','southeast'); end
formatfig
grid on
grid minor
xlim([-1 10])
f2.Position = [326,236,674,561];



end

