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

%% Constants
load('/Users/asper101/Box Sync/Matt Lab Stuff/3rd Year/Electromechanical Modeling/CP_vs_index_fit.mat');
load('BEMTdata_coax_final.mat');

%scale power empirically by steady index
Ct_BEMT = interp1(theta,Ct,11);
Cp_BEMT = interp1(theta,Cp,11);
Cp_up_scale = interp1(index_line,Cp_up,PhaseSync.index_avg(1));
Cp_lo_scale = interp1(index_line,Cp_lo,PhaseSync.index_avg(1));

Cp_up_scale = Cp_up_scale/Cp_BEMT;
Cp_lo_scale = Cp_lo_scale/Cp_BEMT;

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
Thrust_servo.time = Thrust_servo.time - sim_step;

%follower responses
Angle_err_follower.time = Angle_err_follower.time - sim_step; Angle_err_follower.data = Angle_err_follower.data*-1; %reversed reference
Act_angle_follower.time = Act_angle_follower.time - sim_step;
Motor_RPM_follower.time = Motor_RPM_follower.time - sim_step;
Rotor_RPM_follower.time = Rotor_RPM_follower.time - sim_step;
Q_total_follower.time = Q_total_follower.time - sim_step;
Q_motor_follower.time = Q_motor_follower.time - sim_step;
Thrust_follower.time = Thrust_follower.time - sim_step;

%% Plotting
clc;close all;

%Response vs Time
f10 = figure('Name','Response_vs_time');

%Position
subplot(3,1,1)
hold on
plot(Act_angle_servo.time,(Act_angle_servo.data)/GR,'k-','linewidth',1.5)
plot(Act_angle_follower.time,(Act_angle_follower.data)/GR,'k-.','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time+.03,PhaseSync.ref_ang1_avg,PhaseSync.ref_ang1_err,colors{1})
plot_areaerrorbar(PhaseSync.time+.03,PhaseSync.ref_ang2_avg-.5,PhaseSync.ref_ang2_err,colors{2})
ylabel('$\Delta\psi$, deg')
legend('Upper','Lower','Upper','Lower','location','northoutside','orientation','horizontal')
xlim([-.1 1])
formatfig
grid on
grid minor

%Speed
subplot(3,1,2)
hold on
plot(Rotor_RPM_servo.time,Rotor_RPM_servo.data,'k-','linewidth',1.5)
plot(Rotor_RPM_follower.time,Rotor_RPM_follower.data,'k-.','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time+.045,PhaseSync.servo_speed_avg,PhaseSync.servo_speed_err,colors{1})
plot_areaerrorbar(PhaseSync.time+.045,PhaseSync.follower_speed_avg,PhaseSync.follower_speed_err,colors{2})
ylabel('$\Omega$, RPM')
% legend('Upper','Lower','Upper','Lower','location','northeast')
xlim([-.1 1])
formatfig
grid on
grid minor

%Torque
subplot(3,1,3)
hold on
plot(Q_total_servo.time,Q_total_servo.data*Cp_up_scale,'k-','linewidth',1.5)
plot(Q_total_follower.time,Q_total_follower.data*Cp_lo_scale,'k-.','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time,PhaseSync.Torque_outer_avg,PhaseSync.Torque_outer_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.Torque_inner_avg,PhaseSync.Torque_inner_err,colors{2})
ylabel('Q, N$\cdot$m')
% legend('Upper','Lower','Upper','Lower','location','southeast')
xlim([-.1 1])
formatfig
grid on
grid minor
xlabel('Time, s')

f10.Position = [326,236,674,561];

%Estimated Response vs Time
f11 = figure('Name','Estimated_Response_vs_time');

%Position
subplot(3,1,1)
hold on
plot(Act_angle_servo.time,(Act_angle_servo.data)/GR,'k-','linewidth',1.5)
plot(Act_angle_follower.time,(Act_angle_follower.data)/GR,'k-.','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time+.03,PhaseSync.ref_ang1_avg,PhaseSync.ref_ang1_err,colors{1})
plot_areaerrorbar(PhaseSync.time+.03,PhaseSync.ref_ang2_avg-.5,PhaseSync.ref_ang2_err,colors{2})
ylabel('$\Delta\psi$, deg')
legend('Upper','Lower','Upper','Lower','location','northoutside','orientation','horizontal')
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
ylabel('$\Omega$, RPM')
% legend('Upper','Lower','Upper','Lower','location','northeast')
xlim([-.1 1])
formatfig
grid on
grid minor

%Torque
subplot(3,1,3)
hold on
plot(Q_total_servo.time,Q_total_servo.data*Cp_up_scale,'k-','linewidth',1.5)
plot(Q_total_follower.time,Q_total_follower.data*Cp_lo_scale,'k-.','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time,PhaseSync.Q2_est_avg,PhaseSync.Q2_est_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.Q1_est_avg,PhaseSync.Q1_est_err,colors{2})
ylabel('Q, N$\cdot$m')
% legend('Upper','Lower','Upper','Lower','location','southeast')
xlim([-.1 1])
formatfig
grid on
grid minor
xlabel('Time, s')
f11.Position = [326,236,674,561];

%Index
f12 = figure('Name','Index_Angle_vs_Time');
hold on
plot(index.time,index.data,'k-','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time+.03,PhaseSync.index_avg,PhaseSync.index_err,colors{1})
ylabel('$\phi$, deg')
xlabel('Time, s')
xlim([-.1 1])
formatfig
grid on
grid minor

%Ct and Cp vs time
f13 = figure('Name','CT and CP versus Time');
subplot(2,1,1)
hold on
yline(Ct_BEMT/(sigma),'k--','Linewidth',1.5)
plot_areaerrorbar(PhaseSync.time,PhaseSync.cts_up_avg,PhaseSync.cts_up_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.cts_lo_avg,PhaseSync.cts_lo_err,colors{2})
plot_areaerrorbar(PhaseSync.time,PhaseSync.cts_tot_avg,PhaseSync.cts_tot_err,colors{3})
hold off
ylabel('$C_T/\sigma$')
legend('2-bladed','Upper','Lower','Total')
formatfig
xlim([-.1 1])

subplot(2,1,2)
hold on
yline(Cp_BEMT/(sigma),'k--','Linewidth',1.5)
plot_areaerrorbar(PhaseSync.time,PhaseSync.cps_up_avg,PhaseSync.cps_up_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.cps_lo_avg,PhaseSync.cps_lo_err,colors{2})
plot_areaerrorbar(PhaseSync.time,PhaseSync.cps_tot_avg,PhaseSync.cps_tot_err,colors{3})
hold off
xlabel('Time, s')
ylabel('$C_P/\sigma$')
legend('2-bladed','Upper','Lower','Total')
formatfig
xlim([-.1 1])

%Index
f14 = figure('Name','Thrust_vs_Time');
hold on
plot(Thrust_servo.time,Thrust_servo.data,'k--','linewidth',1.5)
plot(Thrust_follower.time,Thrust_follower.data,'k-.','linewidth',1.5)
plot_areaerrorbar(PhaseSync.time+.03,PhaseSync.T_outer_avg,PhaseSync.T_outer_err,colors{1})
plot_areaerrorbar(PhaseSync.time+.03,PhaseSync.T_inner_avg,PhaseSync.T_inner_err,colors{2})
ylabel('Thrust, N')
xlabel('Time, s')
legend('Upper','Lower','Upper','Lower','location','southeast')
xlim([-.1 1])
formatfig
grid on
grid minor
%% Saving
save_dir = uigetdir();
saveas(f10,fullfile(save_dir,f10.Name),'jpg');
saveas(f11,fullfile(save_dir,f11.Name),'jpg');
saveas(f12,fullfile(save_dir,f12.Name),'jpg');
saveas(f13,fullfile(save_dir,f13.Name),'jpg');
saveas(f14,fullfile(save_dir,f14.Name),'jpg');