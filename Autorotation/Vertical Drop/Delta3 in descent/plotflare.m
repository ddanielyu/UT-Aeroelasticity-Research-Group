statesnext = states_theoretical;
rho = physics.rho;
R = rotor.radius;
%% Theta0 plot
figure()
subplot(6,1,1)
hold all
plot(time_vec,rotor.theta0*180/pi);
title('Root collective angle [deg] vs time')

%% Velocity plot
subplot(6,1,2)
hold all
plot(time_vec,statesnext(end,:));
title('Descent velocity vs time')
%% Height plot
% figure()
subplot(6,1,3)
hold all
height = cumtrapz(statesnext(end,:).*dt);
plot(time_vec,height);
title('Height vs time')
%% Diff plot
% figure()
subplot(6,1,4)
hold all
diff1 = diff(statesnext(end,:));
diff2 = diff(statesnext(end,:),2);
plot(time_vec(1:end-1),diff1);
plot(time_vec(1:end-2),diff2);
title('Diffs (dV/dt, d^2V/dt^2)vs time')

%% Thrust plot
subplot(6,1,5)
hold all
plot(time_vec,F(1,:).*statesnext(end-1,:).^2*rho*pi*R^4);
title('Rotor Thrust vs time')
%% RPM plot
subplot(6,1,6)
hold all
plot(time_vec,statesnext(end-1,1:end)*60/2/pi)
title('Rotor RPM vs time')
%% Flap Angle plot
% figure()
% plot(time_vec,statesnext(2,:)*180/pi);
% title('Flap angle vs time')
%}