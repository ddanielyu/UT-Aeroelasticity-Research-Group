function [f1,f2,f3,f4,f5] = plotPhaseSync(PhaseSync)
%{
EDITED ON: 03/17/2022
EDITED BY: MATT ASPER

DETAILS: This update plots estimated motor torque with IQ

%}
%This function plots phase sync data from Dual Motor Testing

%% Plot
load('colors.mat')

%Subplots
f1 = figure('Name','Response_vs_time');
subplot(3,1,1)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.ref_ang1_avg,PhaseSync.ref_ang1_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.ref_ang2_avg,PhaseSync.ref_ang2_err,colors{2})
hold off
ylabel('Rotor Angle, deg')
legend('Upper','Lower','location','southeast')
grid on
grid minor
xlim([-.1 1])

subplot(3,1,2)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.servo_speed_avg,PhaseSync.servo_speed_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.follower_speed_avg,PhaseSync.follower_speed_err,colors{2})
hold off
ylabel('$\Omega_{rotor}$, RPM')
legend('Upper','Lower','location','northeast')
grid on
grid minor
xlim([-.1 1])

subplot(3,1,3)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.Torque_outer_avg,PhaseSync.Torque_outer_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.Torque_inner_avg,PhaseSync.Torque_inner_err,colors{2})
hold off
ylabel('$Q_{rotor}$, $N\cdot m$')
xlabel('Time, s')
legend('Upper','Lower','location','southeast')
formatfig
grid on
xlim([-.1 1])
f1.Position = [326,236,674,561];

%Subplots
f2 = figure('Name','Estimated Response_vs_time');
subplot(3,1,1)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.ref_ang1_avg,PhaseSync.ref_ang1_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.ref_ang2_avg,PhaseSync.ref_ang2_err,colors{2})
ylabel('Rotor Angle, deg')
legend('Upper','Lower','location','southeast')
hold off
grid on
grid minor
xlim([-.1 1])

subplot(3,1,2)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.servo_speed_avg,PhaseSync.servo_speed_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.follower_speed_avg,PhaseSync.follower_speed_err,colors{2})
ylabel('$\Omega_{rotor}$, RPM')
legend('Upper','Lower','location','northeast')
hold off
grid on
grid minor
xlim([-.1 1])

subplot(3,1,3)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.Q2_est_avg,PhaseSync.Q2_est_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.Q1_est_avg,PhaseSync.Q1_est_err,colors{2})
hold off
ylabel('$Q_{rotor}$, $N\cdot m$')
xlabel('Time, s')
legend('Upper','Lower','location','southeast')
formatfig
grid on
xlim([-.1 1])
f2.Position = [326,236,674,561];

%index vs time
f3 = figure('Name','Index Angle versus Time');
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.index_avg,PhaseSync.index_err,colors{1})
hold off
xlabel('Time, s')
ylabel('Rotor Index, $\phi$ [$^\circ$]')
formatfig
xlim([-.1 1])

%Thrust vs time
f4 = figure('Name','Thrust versus Time');
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.T_outer_avg,PhaseSync.T_outer_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.T_inner_avg,PhaseSync.T_inner_err,colors{2})
hold off
xlabel('Time, s')
ylabel('Thrust, N')
legend('Upper','Lower')
formatfig
xlim([-.1 1])

%Ct and Cp vs time
f5 = figure('Name','CT and CP versus Time');
subplot(2,1,1)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.cts_up_avg,PhaseSync.cts_up_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.cts_lo_avg,PhaseSync.cts_lo_err,colors{2})
plot_areaerrorbar(PhaseSync.time,PhaseSync.cts_tot_avg,PhaseSync.cts_tot_err,colors{3})
hold off
ylabel('$C_t/\sigma$')
legend('Upper','Lower','Total')
formatfig
xlim([-.1 1])

subplot(2,1,2)
hold on
plot_areaerrorbar(PhaseSync.time,PhaseSync.cps_up_avg,PhaseSync.cps_up_err,colors{1})
plot_areaerrorbar(PhaseSync.time,PhaseSync.cps_lo_avg,PhaseSync.cps_lo_err,colors{2})
plot_areaerrorbar(PhaseSync.time,PhaseSync.cps_tot_avg,PhaseSync.cps_tot_err,colors{3})
hold off
xlabel('Time, s')
ylabel('$C_p/\sigma$')
legend('Upper','Lower','Total')
formatfig
xlim([-.1 1])

end

