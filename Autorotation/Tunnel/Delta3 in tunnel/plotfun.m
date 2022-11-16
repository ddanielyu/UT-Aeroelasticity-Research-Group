fontname = 'Times';
Vinf = control.Vd;
R = rotor.radius;
rho = physics.rho;
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);
fontsize = 14;
set(0,'defaultaxesfontsize',fontsize);
set(0,'defaulttextfontsize',fontsize);
set(0,'defaulttextinterpreter','tex')

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',2);   % set the default line width to lw
set(0,'defaultLineMarkerSize',12); % set the default line marker size to msz
set(0,'defaultLineLineWidth',2);   % set the default line width to lw
set(0,'defaultErrorBarMarkerSize',8); % set the default line marker size to msz
%% Height plot
figure()
% plot(t,cumtrapz(states_theoretical(end,:))*dt);
% title('Height vs time')
%% Thrust Plot
t = cumtrapz(dpsi(1:end)./states_theoretical(end,1:end));
subplot(3,2,1)
% t0 = 0.412*(Vinf==-5&&tw==1)+0.254*(Vinf==-7&&tw==0)+0.3423*(Vinf==-5&&tw==0);
plot(t,F(1,:).*states_theoretical(end,:).^2*rho*pi*R^4);
title('Rotor Thrust [N] vs time')
% figure()
% plot(cumtrapz(dpsi(1:end)./states_theoretical(end,1:end)),states_theoretical(2,1:end)*180/pi)
% title('Flap angle [\circ] vs time')
subplot(3,2,3)
plot(t,states_theoretical(end,1:end)*60/2/pi)
title('RPM vs time')
% figure()
vh = sqrt(F(1,1:end)*rho*pi*R^4.*states_theoretical(end-1,1:end).^2/2/rho/pi/R^2);
%plot(Vinf./vh,Vinf/2./vh - sqrt(Vinf^2/4./vh.^2-1))
% figure()
% plot(t,stak(2,1:end)*180/pi)
% title('Aerodynamic and inertial flap moments [N] vs time')
% figure();hold all; 
% plot(t,states_theoretical(2,1:end)/max(states_theoretical(2,1:end)));
% plot(t,stak(2,1:end)/max(stak(2,1:end)))
%% Side Forces plot
% figure()
% hold all
% plot(t,movmean(F(end-1,:),10));
% plot(t,movmean(F(end-2,:),10));
% yyaxis right
% plot(t,movmean(atand(F(end-2,:)./F(end-1,:)),10));
% 
% % plot(t,F(end-1,:));
% % plot(t,F(end-2,:));
% title('Side Force vs time')
%% Beta plot
subplot(3,2,5)
hold all
h1 = plot(t,states_theoretical(2,:)*180/pi);
% h2 = plot(t,states_theoretical(4,:)*180/pi);
title('Flap Angles [\circ] vs time')
% xlabel('Time [s]')
% ylabel('Magnitudes of flap and pitch angles [\circ]')
% legend([h1,h2],'Flap angle 1, \beta_1','Flap angle 2, \beta_2','Location','NorthEast')
ylim([0,45])
% xlim([0,3])
box on
%% Phase plot
% figure()
% hold all
% h1 = plot(psi_vec/2/pi,(control.theta1cbar+control.theta1c).*cos(psi_vec)*180/2/pi);
% h2 = plot(psi_vec/2/pi,(states_theoretical(2,:)-states_theoretical(4,:))*180/2/pi);
% h3 = plot(psi_vec/2/pi,(states_theoretical(4,:)-states_theoretical(2,:))*180/2/pi);
% xlabel('Time [s]')
% ylabel('Phase between flap and pitch angles')
% legend([h1,h2,h3],{'Longitudinal Pitch angle, \theta_{1c}'...
%     ,'Flap angle 1,... \beta_1'...
%     ,'Flap angle 2, \beta_2'},'Location','NorthEast')
% ylim([-10,10])
% % xlim([0,3])
% box on
%% BET plots
subplot(3,2,2)
plot(stak(end).Fy1)
title('Sectional Driving Force vs r/R')
subplot(3,2,4)
plot(stak(end).Fz1)
title('Sectional Lifting Force vs r/R')
subplot(3,2,6)
plot(stak(end).alpha*180/pi)
title('Sectional AoA [\circ] vs r/R')