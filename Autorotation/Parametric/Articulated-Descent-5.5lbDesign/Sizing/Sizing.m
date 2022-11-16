CD_vec = 0.9:0.1:1.4;
R_ft = 1.5:.25:3;
R_vec = convlength(R_ft,'ft','m');
S_vec = pi*R_vec.^2;
Weight = convmass(5.5,'lbm','kg')*10;
v = sqrt(Weight./(0.5*CD_vec'*S_vec*1.225));
plot(R_vec,v,'linewidth',2)
xlabel('Rotor Radius [m]')
ylabel('Terminal Descent Velocity')
text(0.35*ones(length(v(:,1)),1),v(:,1),{'C_D = 0.9','C_D = 1.0','C_D = 1.1',...
    'C_D = 1.2','C_D = 1.3','C_D = 1.4'},'FontSize',12)
xlim([0.3,1])
ax = gca;
%% CHANGE PROPERTIES HERE
ax.FontSize = 16;%default font size x for axes tick labels
ax.LabelFontSizeMultiplier = 1.1;%default font size 1.1x for axes labels
ax.FontName = 'Times';%changes font to Times
ax.FontSizeMode = 'manual';%font-size fixed irrespective of axes size
ax.XGrid = 'on';
ax.YGrid = 'on';