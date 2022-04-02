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
dir = uigetdir('/Users/asper101/Box Sync/For Matt/3rd Year/Tripod Stand'); %filename
load(fullfile(dir,'Data.mat'));
GR = input('Gear Ratio: ');

%% Reponses versus time
[sim_file,sim_path] = uigetfile(strcat('/Users/asper101/Box Sync/For Matt/3rd Year/Electromechanical Modeling/'));
load(fullfile(sim_path,sim_file));

sim_step = input('Step Time: ');

index.time = index.time - sim_step;

%servo responses
Angle_err_servo.time = Angle_err_servo.time - sim_step; Angle_err_servo.data = Angle_err_servo.data*-1; %reversed reference
Act_angle_servo.time = Act_angle_servo.time - sim_step;
Motor_RPM_servo.time = Motor_RPM_servo.time - sim_step;
Rotor_RPM_servo.time = Rotor_RPM_servo.time - sim_step;
Q_total_servo.time = Q_total_servo.time - sim_step;
Q_motor_servo.time = Q_motor_servo.time - sim_step;

%follower responses
Angle_err_follower.time = Angle_err_follower.time - sim_step; Angle_err_follower.data = Angle_err_follower.data*-1; %reversed reference
Act_angle_follower.time = Act_angle_follower.time - sim_step;
Motor_RPM_follower.time = Motor_RPM_follower.time - sim_step;
Rotor_RPM_follower.time = Rotor_RPM_follower.time - sim_step;
Q_total_follower.time = Q_total_follower.time - sim_step;
Q_motor_follower.time = Q_motor_follower.time - sim_step;

%% Plotting
clc;close all;

%Response vs Time
f10 = figure('Name','Estimated_Response_vs_time');

%Position
subplot(3,1,1)
hold on
plot(Act_angle_servo.time,(Act_angle_servo.data)/GR,'k-','linewidth',1.5)
plot(Act_angle_follower.time,(Act_angle_follower.data-6)/GR,'k-.','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time+.03,PhaseSync.ref_ang1_avg,PhaseSync.ref_ang1_err,colors{1})
plot_areaerrorbar(PhaseSync.time+.03,PhaseSync.ref_ang2_avg+.4,PhaseSync.ref_ang2_err,colors{2})
ylabel('Rotor Angle, deg')
legend('Upper','Lower','Upper','Lower','location','southeast')
xlim([-.1 1])
formatfig
grid on
grid minor

%Speed
subplot(3,1,2)
hold on
plot(Rotor_RPM_servo.time,Rotor_RPM_servo.data,'k-','linewidth',1.5)
plot(Rotor_RPM_follower.time,Rotor_RPM_follower.data,'k-.','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time+.05,PhaseSync.servo_speed_avg,PhaseSync.servo_speed_err,colors{1})
plot_areaerrorbar(PhaseSync.time+.05,PhaseSync.follower_speed_avg,PhaseSync.follower_speed_err,colors{2})
ylabel('$\Omega_{rotor}$, RPM')
legend('Upper','Lower','Upper','Lower','location','northeast')
xlim([-.1 1])
formatfig
grid on
grid minor

%Torque
subplot(3,1,3)
hold on
plot(Q_total_servo.time,Q_total_servo.data,'k-','linewidth',1.5)
plot(Q_total_follower.time,Q_total_follower.data/1.28,'k-.','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time+.03,PhaseSync.Q1_est_avg,PhaseSync.Q1_est_err,colors{1})
plot_areaerrorbar(PhaseSync.time+.03,PhaseSync.Q2_est_avg,PhaseSync.Q2_est_err,colors{2})
ylabel('$Q_{rotor}$, $N\cdot m$')
legend('Upper','Lower','Upper','Lower','location','southeast')
xlim([-.1 1])
formatfig
grid on
grid minor

%Index
f11 = figure('Name','Index_Angle_vs_Time');
hold on
plot(index.time,index.data,'k-','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time+.03,PhaseSync.index_avg,PhaseSync.index_err,colors{1})
ylabel('Rotor Index, $\phi$ [deg]')
xlabel('Time, s')
xlim([-.1 1])
formatfig
grid on
grid minor

% sgtitle('$+5^\circ$ Angle Offset','Fontsize',24)
f10.Position = [326,236,674,561];