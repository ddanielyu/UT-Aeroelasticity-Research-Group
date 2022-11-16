statesnext = states_theoretical;
rho = physics.rho;
R = rotor.radius;



%% Thrust plot
figure()
subplot(4,1,1)
plot(time_vec,F(1,:).*statesnext(end-1,:).^2*rho*pi*R^4);
ylabel({'T';'[N]'})
format_figure;
% title('Rotor Thrust vs time')
%% RPM plot
subplot(4,1,2)
plot(time_vec,statesnext(end-1,1:end)*60/2/pi)
ylabel({'\Omega','[RPM]'})
format_figure;
% title('Rotor RPM vs time')
%% Velocity plot
subplot(4,1,3)
hold all
plot(time_vec,statesnext(end,:));
ylabel({'V_\infty', '[m/s]'})
format_figure;
% title('Descent velocity vs time')
%% Height plot
% figure()
subplot(4,1,4)
hold all
height = cumtrapz(statesnext(end,:).*dt);
plot(.4+time_vec,convlength(height,'m','ft'));
xlabel('Time [s]')
ylabel({'H', '[ft]'})
format_figure;
% title('Height vs time')

%% Diff plot
% % figure()
% subplot(3,1,3)
% hold all
% diff1 = diff(statesnext(end,:));
% diff2 = diff(statesnext(end,:),2);
% plot(time_vec(1:end-1),diff1);
% plot(time_vec(1:end-2),diff2);
% title('Diffs vs time')
%% Flap Angle plot
figure()
plot(time_vec,statesnext(2,:)*180/pi);
title('Flap angle vs time')
%}

function format_figure()
box on
grid on
ax = gca;
ax.FontSize = 16;%default font size x for axes tick labels
ax.LabelFontSizeMultiplier = 1.1;%default font size 1.1x for axes labels
ax.FontName = 'Times';%changes font to Times
ax.FontSizeMode = 'manual';%font-size fixed irrespective of axes size
end