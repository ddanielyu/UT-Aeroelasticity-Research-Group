% PLOT AC METRICS FROM FILE
% clear; clc; close all; 
% load('zc15_acmetrics.mat')

% COMPILE DATA
for i = 1:length(testdata)
    oaspl(i,1) = testdata{i}(3).oaspl;
    oasplA(i,1) = testdata{i}(3).oasplA;
    oaspl3(i,1) = testdata{i}(3).oaspl3;
    PL(i,1) = testdata{i}(3).PL;
    PNL(i,1) = testdata{i}(3).PNL;
    oasplA3(i,1) = testdata{i}(3).oasplA3;
end

% AVERAGE DATA
uphi = unique(phi);
for i = 1:length(uphi)
    uoaspl(i) = mean(oaspl(uphi(i)==phi));
    uoasplA(i) = mean(oasplA(uphi(i)==phi));
    uoaspl3(i) = mean(oaspl3(uphi(i)==phi));
    uPL(i) = mean(PL(uphi(i)==phi));
    uPNL(i) = mean(PNL(uphi(i)==phi));
    uoasplA3(i) = mean(oasplA3(uphi(i)==phi));
end

% PLOT DATA
[uphi, loc] = sort(uphi);
uoaspl = uoaspl(loc);
uoaspl3 = uoaspl3(loc);
uoasplA = uoasplA(loc);
uoasplA3 = uoasplA3(loc);
uPL = uPL(loc);
uPNL = uPNL(loc);

figure(1)
hold on
plot(uphi, uoaspl-uoaspl(end-1),'o-')
% plot(uphi, uoaspl3-uoaspl3(end-1),'o-')
plot(uphi, uoasplA-uoasplA(end-1),'o-')
% plot(uphi, uoasplA3-uoasplA3(end-1),'o-')
% plot(uphi, uPL-uPL(end-1),'o-')
plot(uphi, uPNL-uPNL(end-1),'o-')

grid on
% legend('OASPL, dB','1/3 OASPL, dB','OASPL, dBA','1/3 OASPL, dBA','Percieved Loudness, dBPL','Percieved Noise Level, dBPNL')
legend('OASPL, dB','OASPL, dBA','Percieved Noise Level, dBPNL')