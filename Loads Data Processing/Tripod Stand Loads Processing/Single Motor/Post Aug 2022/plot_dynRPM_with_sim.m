clc;close all; clear;
load('colors.mat')
%{
This script plots simulink electromechanical model predictions 
on top of experimental data.

1) Load responses (versus time, revolution, thrust, etc.)
2) Run the appropriate section below
3) Select the desired simulink prediction corresponding to the file notes
.xlsx details

Written By: Matt Asper
Date: 2022 Sept 20
%}
f1 = uigetdir('/Users/asper101/Box Sync/Matt Lab Stuff/4th Year/Test Data'); %filename
load(fullfile(f1,'Data.mat'));
GR = input('Gear Ratio: ');

%% Reponses versus time
[sim_file,sim_path] = uigetfile(strcat('/Users/asper101/Box Sync/Matt Lab Stuff/4th Year/Modeling/Electromechanical Model'));
load(fullfile(sim_path,sim_file));

sim_step = 33;
Rotor_RPM.time = Rotor_RPM.time - sim_step;
Q_total.time = Q_total.time - sim_step;
Thrust.time = Thrust.time - sim_step;
nRev.time = nRev.time - sim_step;
Speed_cmd.time = Speed_cmd.time - sim_step;

[~,angle_step] = min(abs(nRev.time));
nRev.data = nRev.data - nRev.data(angle_step);
%% Plotting
clc;close all;

%Rotor Speed
f1 = figure('Name','Rotor Speed');

%versus time
subplot(2,1,1)
hold on
plot(Speed_cmd.time,Speed_cmd.data*7200/GR,'k-','linewidth',1.5)
plot(Rotor_RPM.time,Rotor_RPM.data,'k--','linewidth',1.5)
plot_areaerrorbar(DynRPM.time-.1,DynRPM.Speed_avg/GR,DynRPM.Speed_err/GR,colors{1})
ylabel('$\Omega$, RPM')
xlabel('Time, s')
legend('Commanded','Prediction','Experiment','location','southeast')
xlim([-1 4])
formatfig
grid on
grid minor

%versus revolution
subplot(2,1,2)
hold on
plot(nRev.data/GR,Speed_cmd.data(1:length(nRev.data))*7200/GR,'k-','linewidth',1.5)
plot(nRev.data/GR,Rotor_RPM.data(1:length(nRev.data)),'k--','linewidth',1.5)
plot_areaerrorbar(DynRPM.rev/360/GR-1,DynRPM.Speed_avg/GR,DynRPM.Speed_err/GR,colors{1})
ylabel('$\Omega$, RPM')
xlabel('Rotor Revolution')
xlim([-15 57])
formatfig
grid on
grid minor

f1.Position = [326,236,674,561];


%Torque
f2 = figure('Name','Rotor Torque');

%versus time
subplot(2,1,1)
hold on
plot(Q_total.time,Q_total.data,'k-','linewidth',1.5)
% plot_areaerrorbar(DynRPM.time,DynRPM.Torque_avg*GR,DynRPM.Torque_err*GR,colors{1})
plot_areaerrorbar(DynRPM.time,DynRPM.Q_est_avg*GR,DynRPM.Q_est_err*GR,colors{2})
ylabel('$Q_r$, N$\cdot$m')
legend('Prediction','Experiment','location','southeast')
xlim([-1 4])
formatfig
grid on
grid minor
xlabel('Time, s')
hold off

%versus revolution
subplot(2,1,2)
hold on
plot(nRev.data/GR,Q_total.data(1:length(nRev.data)),'k-','linewidth',1.5)
% plot_areaerrorbar(DynRPM.rev/360/GR,DynRPM.Torque_avg*GR,DynRPM.Torque_err*GR,colors{1})
plot_areaerrorbar(DynRPM.rev/360/GR,DynRPM.Q_est_avg*GR,DynRPM.Q_est_err*GR,colors{2})
ylabel('$Q_r$, $N\cdot m$')
xlim([-15 57])
formatfig
grid on
grid minor
xlabel('Rotor Revolution')
hold off

f2.Position = [326,236,674,561];


%% Saving
save_dir = uigetdir();
saveas(f1,fullfile(save_dir,f1.Name),'jpeg');
saveas(f2,fullfile(save_dir,f2.Name),'jpeg');