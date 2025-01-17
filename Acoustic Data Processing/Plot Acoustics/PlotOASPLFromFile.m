% PLOT OASPL FROM FILE
clear; clc; close all;
addpath('/Users/chloe/Box/Chloe Lab Stuff/2020 Fall Stacked Rotor/Results');

load('acousticdata.mat')
load('testmatrix.mat');
% load('zc15_ac.mat')
% load('zc15_ac_testmat.mat')
load('colors.mat')

%% GET DATA
diffs = testmat.diff;
phis = testmat.phi; 
colls = testmat.coll;
rpms = testmat.rpm;

%% PLOT REFERENCE
loc = contains(testmat.name, 'ref');
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(10)
l = plot(oaspl,'o');
l.MarkerFaceColor = l.Color;
grid on
xlabel('Reference Number')
ylabel('OASPL, dB')
yl = ylim;
hold on
a=plot([0,0],yl,'--');
b=plot([5,5],yl,'--');
c=plot([8,8],yl,'--');
legend([a,b,c],'201117','201118','201119')

%% PLOT ISOLATED ROTOR
% --------------------990--------------------
clear oaspl oasplA
rpm_des = 990; 
phi_des = 2; 
diff_des = 0; 

loc = (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(1)
l = plot(colls(loc),oaspl,'o');
l.MarkerFaceColor = l.Color;
grid on
xlabel('Collective')
ylabel('OASPL, dB')
title('Isolated Rotor')

%% --------------------1200--------------------
clear oaspl oasplA
rpm_des = 1200; 
phi_des = 2; 
diff_des = 0; 

loc = (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(1)
hold on
l = plot(colls(loc),oaspl,'o');
l.MarkerFaceColor = l.Color;
legend('990 RPM','1200 RPM','location','northwest')

% --------------------1200 201118--------------------
clear oaspl oasplA
rpm_des = 1200; 
phi_des = 2; 
diff_des = 0; 
date_des = '201118';

loc = contains(testmat.name, '201118');
loc = loc & (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(5)
[a,b] = sort(colls(loc));
l = plot(a,oaspl(b),'o-');
l.MarkerFaceColor = l.Color;
grid on
xlabel('Collective')
ylabel('dB')
title('1200 RPM')

%% PLOT SINGLE ROTOR
% --------------------990--------------------
clear oaspl oasplA
rpm_des = 990; 
phi_des = 90; 
diff_des = 0; 

loc = (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(2)
l = plot(colls(loc),oaspl,'o');
l.MarkerFaceColor = l.Color;
grid on
xlabel('Collective')
ylabel('OASPL, dB')
title('4-Bladed Rotor')

% --------------------1200--------------------
clear oaspl oasplA
rpm_des = 1200; 
phi_des = 90; 
diff_des = 0; 

loc = (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(2)
hold on
l = plot(colls(loc),oaspl,'o');
l.MarkerFaceColor = l.Color;
legend('990 RPM','1200 RPM','location','northwest')

% --------------------1200 201118--------------------
clear oaspl oasplA
rpm_des = 1200; 
phi_des = 90; 
diff_des = 0; 
date_des = '201118';

loc = contains(testmat.name, '201118');
loc = loc & (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(5)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oaspl(b),'o-');
l.MarkerFaceColor = l.Color;
legend('2-Bladed','4-Bladed','location','northwest')

% --------------------1200 Differential--------------------
clear oaspl oasplA
rpm_des = 1200; 
phi_des = 90; 
coll_des = 12; 
date_des = '201118';

loc = contains(testmat.name, '201118');
loc = loc & (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (colls == coll_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(3)
hold on
[a,b] = sort(diffs(loc));
l = plot(a,oaspl(b),'o-');
l.MarkerFaceColor = l.Color;
grid on
xlabel('Differential')
ylabel('OASPL, dB')
title('4-Bladed Rotor')

%% PLOT PHIS
% --------------------90--------------------
clear oaspl oasplA
rpm_des = 1200; 
phi_des = 90; 
diff_des = 0; 
% date_des = '201118';

% loc = contains(testmat.name, '201118');
loc = (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(4)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oaspl(b),'o');
l.MarkerFaceColor = l.Color;
grid on
xlabel('Collective')
ylabel('OASPL, dB')

figure(6)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oasplA(b),'o');
l.MarkerFaceColor = l.Color;
grid on
xlabel('Collective')
ylabel('OASPLA, dB')
% --------------------28.125--------------------
clear oaspl oasplA
rpm_des = 1200; 
phi_des = 28.125; 
diff_des = 0; 

loc = (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(4)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oaspl(b),'o');
l.MarkerFaceColor = l.Color;

figure(6)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oasplA(b),'o');
l.MarkerFaceColor = l.Color;
% --------------------39.375--------------------
clear oaspl oasplA
rpm_des = 1200; 
phi_des = 39.375; 
diff_des = 0; 

loc = (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(4)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oaspl(b),'o');
l.MarkerFaceColor = l.Color;

figure(6)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oasplA(b),'o');
l.MarkerFaceColor = l.Color;
% --------------------45--------------------
clear oaspl oasplA
rpm_des = 1200; 
phi_des = 45; 
diff_des = 0; 

loc = (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(4)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oaspl(b),'o');
l.MarkerFaceColor = l.Color;

figure(6)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oasplA(b),'o');
l.MarkerFaceColor = l.Color;
% --------------------50.625--------------------
clear oaspl oasplA
rpm_des = 1200; 
phi_des = 50.625; 
diff_des = 0; 

loc = (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(4)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oaspl(b),'o');
l.MarkerFaceColor = l.Color;

figure(6)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oasplA(b),'o');
l.MarkerFaceColor = l.Color;
% --------------------67.5--------------------
clear oaspl oasplA
rpm_des = 1200; 
phi_des = 67.5; 
diff_des = 0; 

loc = (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (phis == phi_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
end
figure(4)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oaspl(b),'o');
l.MarkerFaceColor = l.Color;
l = legend('90','28.125','39.375','45','50.625','67.5','location','northwest');
title(l,'\phi')

figure(6)
hold on
[a,b] = sort(colls(loc));
l = plot(a,oasplA(b),'o');
l.MarkerFaceColor = l.Color;

l = legend('90','28.125','39.375','45','50.625','67.5','location','northwest');
title(l,'\phi')

%% PLOT VS PHI
clear oaspl oasplA
rpm_des = 1200; 
coll_des = 10; 
diff_des = 0; 

loc = (contains(testmat.name, '201118'))|(contains(testmat.name, '201119'));
loc = (phis ~= 2) & loc & (rpms > rpm_des*.98) & (rpms < rpm_des*1.02) & (colls == coll_des) & (diffs == diff_des);
plotdata = {data{loc}};
for i = 1:length(plotdata)
    oaspl(i,1) = plotdata{i}(3).oaspl;
    oasplA(i,1) = plotdata{i}(3).oasplA;
    tl(i,1) = plotdata{i}(3).dBtl;
    bb(i,1) = plotdata{i}(3).dBbb;    
end

phi_plot = phis(loc);
phi_uni = unique(phi_plot);
for i = 1:length(phi_uni)
    loc = phi_plot==phi_uni(i);
    oaspls(i,1) = mean(oaspl(loc));
    oasplAs(i,1) =mean(oasplA(loc));
    bbs(i,1) = mean(bb(loc));
    tls(i,1) = mean(tl(loc));
end

figure(7)
    hold on
    [a,b] = sort(phi_uni);
    l = plot(a,oaspls(b),'ko-');
    l.MarkerFaceColor = l.Color;
    
    hold on
    [a,b] = sort(phi_uni);
    l = plot(a,tls(b),'o-','color',colors{1});
    l.MarkerFaceColor = l.Color;
    
    hold on
    [a,b] = sort(phi_uni);
    l = plot(a,bbs(b),'o-','color',colors{2});
    l.MarkerFaceColor = l.Color;

    xlabel('Azimuthal Spacing, \phi')
    ylabel('SPL, dB')
    legend('Exp Total','Exp T+L','Exp Broadband')

figure(8)
    hold on
    [a,b] = sort(phis(loc));
    l = plot(a,oasplA(b),'ko-');
    l.MarkerFaceColor = l.Color;
    xlabel('Azimuthal Spacing, \phi')
    ylabel('OASPLA, db')

