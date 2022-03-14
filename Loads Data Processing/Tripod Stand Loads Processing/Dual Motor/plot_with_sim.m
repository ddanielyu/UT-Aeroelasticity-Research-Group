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
f1 = uigetfile(); %filename
f1 = openfig(f1);

%% Reponses versus time
[sim_file,sim_path] = uigetfile(strcat('/Users/asper101/Box Sync/For Matt/3rd Year/Electromechanical Modeling/'));
load(fullfile(sim_path,sim_file));

sim_step = 10;
Angle_err.time = Angle_err.time - sim_step; Angle_err.data = Angle_err.data*-1; %reversed reference
Act_angle.time = Act_angle.time - sim_step;
Motor_RPM.time = Motor_RPM.time - sim_step;
Q_total.time = Q_total.time - sim_step;

figure(f6)
hold on
subplot(3,1,1)
plot(Angle_err.time,Angle_err.data,'k-','linewidth',1.5)

subplot(3,1,2)
plot(Motor_RPM.time,Motor_RPM.data,'k-','linewidth',1.5)

subplot(3,1,3)
plot(Q_total.time,Q_total.data/1.2,'k-','linewidth',1.5)
