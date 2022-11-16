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
% figure()
% plot(t,cumtrapz(states_theoretical(end,:))*dt);
% title('Height vs time')
%% Thrust Plot
t = cumtrapz(dpsi(1:end)./states_theoretical(end,1:end));
figure()
% t0 = 0.412*(Vinf==-5&&tw==1)+0.254*(Vinf==-7&&tw==0)+0.3423*(Vinf==-5&&tw==0);
plot(t,convmass(F(1,:).*states_theoretical(end,:).^2*rho*pi*R^4/9.8,'kg','lbm'));
title('Rotor Thrust [lb] vs time')
% figure()
% plot(cumtrapz(dpsi(1:end)./states_theoretical(end,1:end)),states_theoretical(2,1:end)*180/pi)
% title('Flap angle [\circ] vs time')
figure()
plot(t,states_theoretical(end,1:end)*60/2/pi)
title('RPM vs time')

%% Thrust Plot 2
t = cumtrapz(dpsi(1:end)./states_theoretical(end,1:end));
figure()
% t0 = 0.412*(Vinf==-5&&tw==1)+0.254*(Vinf==-7&&tw==0)+0.3423*(Vinf==-5&&tw==0);
[~,init] = min(abs(t-1.955));
plot(t(init:end)-t(init),...
    convmass(F(1,init:end).*states_theoretical(end,init:end).^2*rho*pi*R^4/9.8,'kg','lbm'),'k');
title('Rotor Thrust [lb] vs time')

figure()
% t0 = 0.412*(Vinf==-5&&tw==1)+0.254*(Vinf==-7&&tw==0)+0.3423*(Vinf==-5&&tw==0);
plot(t(init:end)-t(init),(control.theta0bar+control.theta0(init:end))*180/pi,'k')
title('Collective Input vs Time')
xlabel('Time [s]')
ylabel('Collective Input, \theta_0 [\circ]')