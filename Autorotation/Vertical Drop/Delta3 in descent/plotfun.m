statesnext = states_theoretical;
rho = physics.rho;
R = rotor.radius;

%% Velocity plot
figure()
plot(time_vec,statesnext(end,:));
title('Descent velocity vs time')
%% Height plot
figure()
height = cumtrapz(statesnext(end,:).*dt);
plot(time_vec,height);
title('Height vs time')
%% Thrust plot
figure()
plot(time_vec,F(1,:).*statesnext(end-1,:).^2*rho*pi*R^4);
title('Rotor Thrust vs time')
%% RPM plot
figure()
plot(time_vec,statesnext(end-1,1:end)*60/2/pi)
title('Rotor RPM vs time')
%% Flap Angle plot
figure()
plot(time_vec,statesnext(2,:)*180/pi);
title('Flap angle vs time')
