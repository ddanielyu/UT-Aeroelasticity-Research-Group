RPM = 840; 
k = 5;

prev = RPM/60; 

dt = 1/10000;
tvec = dt:dt:length(StreamData.ax{k})*dt;
[fvec, magVsx, psdVs, dftVs2] = ffind_dft(tvec, StreamData.ax{k});
loc = (fvec>prev*0.8) & (fvec<prev*1.2); 
ax = max(magVsx(loc)); 

[fvec, magVsy, psdVs, dftVs2] = ffind_dft(tvec, StreamData.ay{k});
ay = max(magVsy(loc)); 

figure(1)
subplot(2,1,1)
hold on
plot(fvec, magVsx,'.-')
ylabel('a_x [g]')
ylim([0,0.1])
xlim([0,100])
subplot(2,1,2)
hold on
plot(fvec, magVsy,'.-')
ylim([0,0.1])
xlim([0,100])
ylabel('a_y [g]')
hold off
xlabel('Freq [Hz]')

a = sqrt(ax.^2 + ay.^2)

figure(2)
subplot(2,1,1)
hold on
plot(RevData.avg_ax{k},'.-')
ylabel('a_x [g]')
subplot(2,1,2)
hold on
plot(RevData.avg_ay{k},'.-')
ylabel('a_y [g]')
hold off
xlabel('azimuth')


% bandpass(RevData.avg_ax{k},[prev*0.8,prev*1.2],10000)