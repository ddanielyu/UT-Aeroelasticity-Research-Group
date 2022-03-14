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
Date: 2022 Feb 24
%}
f1 = uigetdir('/Users/asper101/Box Sync/For Matt/3rd Year/Tripod Stand'); %filename
load(fullfile(f1,'Data.mat'));
GR = input('Gear Ratio: ');

%% Reponses versus time
[sim_file,sim_path] = uigetfile(strcat('/Users/asper101/Box Sync/For Matt/3rd Year/Electromechanical Modeling/'));
load(fullfile(sim_path,sim_file));

sim_step = 10;
Angle_err.time = Angle_err.time - sim_step; Angle_err.data = Angle_err.data*-1; %reversed reference
Act_angle.time = Act_angle.time - sim_step;
Motor_RPM.time = Motor_RPM.time - sim_step;
Q_total.time = Q_total.time - sim_step;
Q_motor.time = Q_motor.time - sim_step;

%% Plotting
clc;close all;

%Response vs Time
f1 = figure('Name','Response_vs_time');

%Position
subplot(3,1,1)
hold on
plot(Act_angle.time,Act_angle.data,'k-','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time{end}-.02,PhaseSync.ref_ang_avg{end},PhaseSync.ref_ang_err{end},colors{1})
ylabel('Motor Angle, deg')
legend('Prediction','Experiment','location','southeast')
xlim([-.1 1])
formatfig
grid on
grid minor

%Speed
subplot(3,1,2)
hold on
plot(Motor_RPM.time,Motor_RPM.data,'k-','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time{end}-.02,PhaseSync.Speed_avg{end},PhaseSync.Speed_err{end},colors{1})
ylabel('$\Omega_{motor}$, RPM')
legend('Prediction','Experiment','location','northeast')
xlim([-.1 1])
formatfig
grid on
grid minor

%Torque
subplot(3,1,3)
hold on
plot(Q_motor.time,Q_motor.data,'k-','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time{end}-.02,PhaseSync.Q_est_avg{end},PhaseSync.Q_est_err{end},colors{2})
ylabel('$Q_{motor}$, $N\cdot m$')
legend('Prediction','Experiment','location','southeast')
xlim([-.1 1])
formatfig
grid on
grid minor

sgtitle('$+5^\circ$ Angle Offset','Fontsize',24)
f1.Position = [326,236,674,561];