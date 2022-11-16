clearvars;
main_flag = 1;
%% Iterables
theta0_vec = -8:2:10;
theta_tw_vec = -9:3:6;
AR_vec = 6:2:12;
delta3_vec = (20:10:50)*pi/180;

%% Twist Parameteric
AR = .22/0.025;
delta3 = 30*pi/180;
theta_tw = 0;
theta0 = 0;
ind = 0;
for theta_tw = theta_tw_vec
    ind = ind+1
    main
    diff1 = diff(states_theoretical(end,:));
    diff2 = diff(states_theoretical(end,:),2);
    find_ss(ind,:) = find((states_theoretical(end,1:end-1)-1.01*states_theoretical(end,end))./...
        (states_theoretical(end,2:end)-1.01*states_theoretical(end,end))<0,2,'first');
    ind_ss = find_ss(ind,2);
    t_tw_ss(ind) = time_vec(ind_ss);
    v_tw_ss(ind) = states_theoretical(end,ind_ss);
    h_tw_ss(ind) = trapz(states_theoretical(end,1:ind_ss).*dt(1:ind_ss));
end
% Plot
figure();hold all;title('Performance of Decelerator at Different Twist Gradients')
plot(theta_tw_vec,t_tw_ss/max(abs(t_tw_ss)),theta_tw_vec,-v_tw_ss/max(abs(v_tw_ss)),...
    theta_tw_vec,-h_tw_ss/max(abs(h_tw_ss)),'linewidth',2);
xlabel('Twist angle between root and tip [\circ]');
ylabel('Scale with respect to maximum')
legend(strcat('Time to SS [s] scaled by maximum = ', num2str(max(abs(t_tw_ss)))),...
    strcat('Velocity to SS [m/s] scaled by maximum = ', num2str(max(abs(v_tw_ss)))),...
    strcat('Height descended to SS [m] scaled by maximum = ', num2str(max(abs(h_tw_ss)))))


% ands=(diff1(2:end-1)-diff2(2:end))/(diff1(1:end-2)-diff2(1:end-1))


%% Aspect Ratio Parametric
AR = .22/0.025;
delta3 = 30*pi/180;
theta_tw = 0;
theta0 = 0;
ind = 0;
for AR = AR_vec
    ind = ind+1
    main
    diff1 = diff(states_theoretical(end,:));
    diff2 = diff(states_theoretical(end,:),2);
    find_ss(ind,:) = find((states_theoretical(end,1:end-1)-1.01*states_theoretical(end,end))./...
        (states_theoretical(end,2:end)-1.01*states_theoretical(end,end))<0,2,'first');
    %     find((diff1(2:end-1)-diff2(2:end))./...
    %         (diff1(1:end-2)-diff2(1:end-1))<0,2,'first');
    ind_ss = find_ss(ind,2);
    t_AR_ss(ind) = time_vec(ind_ss);
    v_AR_ss(ind) = states_theoretical(end,ind_ss);
    h_AR_ss(ind) = trapz(states_theoretical(end,1:ind_ss).*dt(1:ind_ss));
    
end
% Plot
figure();hold all;title('Performance of Decelerator at Different Aspect Ratios')
plot(AR_vec,t_AR_ss/max(abs(t_AR_ss)),AR_vec,-v_AR_ss/max(abs(v_AR_ss)),...
    AR_vec,-h_AR_ss/max(abs(h_AR_ss)),'linewidth',2);
xlabel('Aspect Rario');
ylabel('Scale with respect to maximum')
legend(strcat('Time to SS [s] scaled by maximum = ', num2str(max(abs(t_AR_ss)))),...
    strcat('Velocity to SS [m/s] scaled by maximum = ', num2str(max(abs(v_AR_ss)))),...
    strcat('Height descended to SS [m] scaled by maximum = ', num2str(max(abs(h_AR_ss)))))

%% Delta3 Parametric
AR = .22/0.025;
delta3 = 30*pi/180;
theta_tw = 0;
theta0 = 0;
ind = 0;
for delta3 = delta3_vec
    ind = ind+1
    main
    diff1 = diff(states_theoretical(end,:));
    diff2 = diff(states_theoretical(end,:),2);
    find_ss(ind,:) = find((states_theoretical(end,1:end-1)-1.01*states_theoretical(end,end))./...
        (states_theoretical(end,2:end)-1.01*states_theoretical(end,end))<0,2,'first');
    ind_ss = find_ss(ind,2);
    t_d3_ss(ind) = time_vec(ind_ss);
    v_d3_ss(ind) = states_theoretical(end,ind_ss);
    h_d3_ss(ind) = trapz(states_theoretical(end,1:ind_ss).*dt(1:ind_ss));
    
end
% Plot
figure();hold all;title('Performance of Decelerator at Different \delta_3 Angles')
plot(delta3_vec*180/pi,t_d3_ss/max(abs(t_d3_ss)),delta3_vec*180/pi,-v_d3_ss/max(abs(v_d3_ss)),...
    delta3_vec*180/pi,-h_d3_ss/max(abs(h_d3_ss)),'linewidth',2);
xlabel('\delta_3 angle [\circ]');
ylabel('Scale with respect to maximum')
legend(strcat('Time to SS [s] scaled by maximum = ', num2str(max(abs(t_d3_ss)))),...
    strcat('Velocity to SS [m/s] scaled by maximum = ', num2str(max(abs(v_d3_ss)))),...
    strcat('Height descended to SS [m] scaled by maximum = ', num2str(max(abs(h_d3_ss)))))

%% Root Pitch Parametric
AR = .22/0.025;
delta3 = 30*pi/180;
theta_tw = 0;
theta0 = 0;
ind = 0;
for theta0 = theta0_vec*pi/180
    ind = ind+1
    main
    diff1 = diff(states_theoretical(end,:));
    diff2 = diff(states_theoretical(end,:),2);
    find_ss(ind,:) = find((states_theoretical(end,1:end-1)-1.01*states_theoretical(end,end))./...
        (states_theoretical(end,2:end)-1.01*states_theoretical(end,end))<0,2,'first');
    ind_ss = find_ss(ind,2);
    t_0_ss(ind) = time_vec(ind_ss);
    v_0_ss(ind) = states_theoretical(end,ind_ss);
    h_0_ss(ind) = trapz(states_theoretical(end,1:ind_ss).*dt(1:ind_ss));
    
end
% Plot
figure();hold all;title('Performance of Decelerator at Different Root Pitch Angles')
plot(theta0_vec,t_0_ss/max(abs(t_0_ss)),theta0_vec,-v_0_ss/max(abs(v_0_ss)),...
    theta0_vec,-h_0_ss/max(abs(h_0_ss)),'linewidth',2);
xlabel('Root pitch angle [\circ]');
ylabel('Scale with respect to maximum')
legend(strcat('Time to SS [s] scaled by maximum = ', num2str(max(abs(t_0_ss)))),...
    strcat('Velocity to SS [m/s] scaled by maximum = ', num2str(max(abs(v_0_ss)))),...
    strcat('Height descended to SS [m] scaled by maximum = ', num2str(max(abs(h_0_ss)))))

%% Carpet Parametric - delta3, theta0
AR = .22/0.025;
delta3 = 30*pi/180;
theta_tw = -5;
theta0 = 0;
ind = 0;
for delta3 = delta3_vec
    for theta0 = theta0_vec*pi/180
        ind = ind+1
        main
        diff1 = diff(states_theoretical(end,:));
        diff2 = diff(states_theoretical(end,:),2);
        find_ss = find((states_theoretical(end,1:end-1)-1.01*states_theoretical(end,end))./...
            (states_theoretical(end,2:end)-1.01*states_theoretical(end,end))<0,2,'first');
        ind_ss = find_ss(2);
        t_d3_ss(ind) = time_vec(ind_ss);
        v_d3_ss(ind) = states_theoretical(end,ind_ss);
        h_d3_ss(ind) = trapz(states_theoretical(end,1:ind_ss).*dt(1:ind_ss));
        
    end
end
% Plot
t_d3_ss = reshape(t_d3_ss,[length(theta0_vec),length(delta3_vec)]);
v_d3_ss = reshape(v_d3_ss,[length(theta0_vec),length(delta3_vec)]);
h_d3_ss = reshape(h_d3_ss,[length(theta0_vec),length(delta3_vec)]);

[Amat, Bmat] = meshgrid(theta0_vec,delta3_vec*180/pi);

figure();hold all;title('Time to Steady State for Different \delta_3 and \theta_0 Angles')
contourf(Amat,Bmat,t_d3_ss')
xlabel('\theta_0 angle [\circ]');
ylabel('\delta_3 angle [\circ]')
colorbar

figure();hold all;title('Steady State Velocity for Different \delta_3 and \theta_0 Angles')
contourf(Amat,Bmat,v_d3_ss')
xlabel('\theta_0 angle [\circ]');
ylabel('\delta_3 angle [\circ]')
colorbar

figure();hold all;title('Height Descended till Steady State for Different \delta_3 and \theta_0 Angles')
contourf(Amat,Bmat,h_d3_ss')
xlabel('\theta_0 angle [\circ]');
ylabel('\delta_3 angle [\circ]')
colorbar

%% Carpet Parametric - delta3, theta_tw
AR = .22/0.025;
delta3 = 30*pi/180;
theta_tw = 0;
theta0 = 5;
ind = 0;
for delta3 = delta3_vec
    for theta_tw = theta_tw_vec
        ind = ind+1
        main
        diff1 = diff(states_theoretical(end,:));
        diff2 = diff(states_theoretical(end,:),2);
        find_ss = find((states_theoretical(end,1:end-1)-1.01*states_theoretical(end,end))./...
            (states_theoretical(end,2:end)-1.01*states_theoretical(end,end))<0,2,'first');
        ind_ss = find_ss(2);
        t_d3_ss(ind) = time_vec(ind_ss);
        v_d3_ss(ind) = states_theoretical(end,ind_ss);
        h_d3_ss(ind) = trapz(states_theoretical(end,1:ind_ss).*dt(1:ind_ss));
        
    end
end
% Plot
t_d3_ss = reshape(t_d3_ss,[length(theta_tw_vec),length(delta3_vec)]);
v_d3_ss = reshape(v_d3_ss,[length(theta_tw_vec),length(delta3_vec)]);
h_d3_ss = reshape(h_d3_ss,[length(theta_tw_vec),length(delta3_vec)]);

[Amat, Bmat] = meshgrid(theta_tw_vec,delta3_vec*180/pi);

figure();hold all;title('Time to Steady State for Different \delta_3 and \theta_{tw} Angles')
contourf(Amat,Bmat,t_d3_ss')
xlabel('\theta_{tw} angle [\circ]');
ylabel('\delta_3 angle [\circ]')
colorbar

figure();hold all;title('Steady State Velocity for Different \delta_3 and \theta_{tw} Angles')
contourf(Amat,Bmat,v_d3_ss')
xlabel('\theta_{tw} angle [\circ]');
ylabel('\delta_3 angle [\circ]')
colorbar

figure();hold all;title('Height Descended till Steady State for Different \delta_3 and \theta_{tw} Angles')
contourf(Amat,Bmat,h_d3_ss')
xlabel('\theta_{tw} angle [\circ]');
ylabel('\delta_3 angle [\circ]')
colorbar


%% Carpet Parametric - theta0, theta_tw
AR = .22/0.025;
delta3 = 30*pi/180;
theta_tw = 0;
theta0 = 0;
ind = 0;
for theta0 = theta0_vec*pi/180
    for theta_tw = theta_tw_vec
        ind = ind+1
        main
        diff1 = diff(states_theoretical(end,:));
        diff2 = diff(states_theoretical(end,:),2);
        find_ss = find((states_theoretical(end,1:end-1)-1.01*states_theoretical(end,end))./...
            (states_theoretical(end,2:end)-1.01*states_theoretical(end,end))<0,2,'first');
        ind_ss = find_ss(2);
        t_d3_ss(ind) = time_vec(ind_ss);
        v_d3_ss(ind) = states_theoretical(end,ind_ss);
        h_d3_ss(ind) = trapz(states_theoretical(end,1:ind_ss).*dt(1:ind_ss));
        
    end
end
% Plot
t_d3_ss = reshape(t_d3_ss,[length(theta_tw_vec),length(theta0_vec)]);
v_d3_ss = reshape(v_d3_ss,[length(theta_tw_vec),length(theta0_vec)]);
h_d3_ss = reshape(h_d3_ss,[length(theta_tw_vec),length(theta0_vec)]);

[Amat, Bmat] = meshgrid(theta_tw_vec,theta0_vec);

figure();hold all;title('Time to Steady State for Different \theta_0 and \theta_{tw} Angles')
contourf(Amat,Bmat,t_d3_ss')
xlabel('\theta_{tw} angle [\circ]');
ylabel('\theta_0 angle [\circ]')
colorbar

figure();hold all;title('Steady State Velocity for Different \theta_0 and \theta_{tw} Angles')
contourf(Amat,Bmat,v_d3_ss')
xlabel('\theta_{tw} angle [\circ]');
ylabel('\theta_0 angle [\circ]')
colorbar

figure();hold all;title('Height Descended till Steady State for Different \theta_0 and \theta_{tw} Angles')
contourf(Amat,Bmat,h_d3_ss')
xlabel('\theta_{tw} angle [\circ]');
ylabel('\theta_0 angle [\circ]')
colorbar


%% Carpet Parametric - theta0, theta_tw, delta_3
AR = .22/0.025;
delta3 = 30*pi/180;
theta_tw = 0;
theta0 = 0;
ind = 0;
for theta0 = theta0_vec*pi/180
    for theta_tw = theta_tw_vec
        for delta_3 = delta3_vec
            ind = ind+1
            main
            %             diff1 = diff(states_theoretical(end,:));
            %             diff2 = diff(states_theoretical(end,:),2);
            find_ss = find((states_theoretical(end,:)<0.95*states_theoretical(end,end)),1,'last');
            %             [~,find_ss] = min(states_theoretical(end,1:end-1));
            %             [~,find_ss] = find((states_theoretical(end,1:end-1)-1.01*states_theoretical(end,end))./...
            %                 (states_theoretical(end,2:end)-1.01*states_theoretical(end,end))<0,2,'first');
            ind_ss = find_ss;%(2);
            t_d3_ss(ind) = time_vec(ind_ss);
            v_d3_ss(ind) = states_theoretical(end,ind_ss);
            h_d3_ss(ind) = trapz(states_theoretical(end,1:ind_ss).*dt(1:ind_ss));
            CD_d3_ss(ind) = F(1,ind_ss)*rotor.radius^2*states_theoretical(end-1,ind_ss)^2/0.5/states_theoretical(end,ind_ss)^2;
        end
    end
end
%%
ind = 0;
for theta0 = theta0_vec*pi/180
    for theta_tw = theta_tw_vec
        for delta_3 = delta3_vec
            ind = ind+1
            theta0_scatter(ind) = theta0;
            theta_tw_scatter(ind) = theta_tw;
            delta3_scatter(ind) = delta_3;
        end
    end
end
%% Plot
t1_d3_ss = reshape(t_d3_ss,[length(theta_tw_vec),length(theta0_vec),length(delta3_vec)]);
v1_d3_ss = reshape(v_d3_ss,[length(theta_tw_vec),length(theta0_vec),length(delta3_vec)]);
h1_d3_ss = reshape(h_d3_ss,[length(theta_tw_vec),length(theta0_vec),length(delta3_vec)]);
CD1_d3_ss = reshape(CD_d3_ss,[length(theta_tw_vec),length(theta0_vec),length(delta3_vec)]);

[Amat, Bmat] = meshgrid(theta_tw_vec,theta0_vec(2:end));%,length(delta3_vec));
%%
figure();
map = 'default';
line_styles = {'-','--','-.',':'};
for i = 1:length(delta3_vec)
    surface(Amat',Bmat',squeeze(CD1_d3_ss(:,2:end,i)),'LineStyle',line_styles{i},'FaceColor','interp')
    colormap(map);
end
caxis([.9 1.1]);
view(0, 90);
xlabel('\theta_{tw} [deg. from root to tip]')
ylabel('\theta_0 [deg.]')
legend('\delta_3 = 20\circ','\delta_3 = 30\circ','\delta_3 = 40\circ','\delta_3 = 50\circ')
clrbr = colorbar;
clrbr.Label.String = 'Time to steady state [s]';

%%
figure();
map = 'jet';
line_styles = {'-','--','-.',':'};
for i = 1:length(delta3_vec)
    surface(Amat',Bmat',-squeeze(v_d3_ss(:,2:end,i)),'LineStyle',line_styles{i})
    colormap(map);
end
% caxis([4.8 5]);
view(0, -90);
xlabel('\theta_{tw} [deg. from root to tip]')
ylabel('\theta_0 [deg.]')
legend('\delta_3 = 20\circ','\delta_3 = 30\circ','\delta_3 = 40\circ','\delta_3 = 50\circ')
clrbr = colorbar;
clrbr.Label.String = 'Steady state descent velocity [m/s]';

%%
figure();
map = 'default';
line_styles = {'-','--','-.',':'};
for i = 1:length(delta3_vec)
    surface(Amat',Bmat',-squeeze(h_d3_ss(:,2:end,i)),'LineStyle',line_styles{i})
    colormap(map);
end
% caxis([13 17]);
view(0, -90);
xlabel('\theta_{tw} [deg. from root to tip]')
ylabel('\theta_0 [deg.]')
legend('\delta_3 = 20\circ','\delta_3 = 30\circ','\delta_3 = 40\circ','\delta_3 = 50\circ')
clrbr = colorbar;
clrbr.Label.String = 'Steady state height descended [m]';

%%
%ax=-5,el=13
figure()
scatter3(theta0_scatter*180/pi,theta_tw_scatter,delta3_scatter*180/pi,CD_d3_ss*50,CD_d3_ss,'filled')
box on
xlabel('Root Collective, \theta_0 [\circ]')
ylabel('Blade Twist, \theta_{tw} [\circ]')
ylim([-6,6])
yticks(-6:3:6)
zlabel('\delta_3 [\circ]')
ax = gca;
ax.BoxStyle = 'full';
ax.FontSize = 16;%default font size x for axes tick labels
ax.LabelFontSizeMultiplier = 1.1;%default font size 1.1x for axes labels
ax.FontName = 'Times';%changes font to Times
ax.FontSizeMode = 'manual';%font-size fixed irrespective of axes size
ax.Colormap = hsv;
c = colorbar;
c.Label.String = 'Decelerator Drag Coefficient, C_D';
grid3 on
view([-5 13])

%%
%ax=-5,el=13
figure()
scatter3(theta0_scatter*180/pi,theta_tw_scatter,delta3_scatter*180/pi,-h_d3_ss,h_d3_ss,'filled')
box on
xlabel('Root Collective, \theta_0 [\circ]')
ylabel('Blade Twist, \theta_{tw} [\circ]')
ylim([-6,6])
yticks(-6:3:6)
zlabel('\delta_3 [\circ]')
ax = gca;
ax.BoxStyle = 'full';
ax.FontSize = 16;%default font size x for axes tick labels
ax.LabelFontSizeMultiplier = 1.1;%default font size 1.1x for axes labels
ax.FontName = 'Times';%changes font to Times
ax.FontSizeMode = 'manual';%font-size fixed irrespective of axes size
ax.Colormap = hsv;
c = colorbar;
c.Label.String = 'Height Descended until Steady State, [s]';
grid3 on
view([-5 13])