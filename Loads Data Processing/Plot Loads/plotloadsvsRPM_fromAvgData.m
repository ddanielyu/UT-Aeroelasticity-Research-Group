clear RPMs Fzo Fzi Mzo Mzi
RPMs = MeanData.RPMs;
for i = 1:length(RPMs)
Fzo(i) = mean([StreamData.Fz_outer{i}]);
Fzi(i) = mean([StreamData.Fy_outer{i}]);
Mzo(i) = -mean([StreamData.Mz_outer{i}]);
Mzi(i) = mean([StreamData.Mz_inner{i}]);
end
figure(6)
hold on
plot(RPMs, Fzo*0.224809,'o')
plot(RPMs, Fzi*0.224809,'o')
xlabel('RPM')
ylabel('Thrust, lbf')
figure(2)
hold on
plot(RPMs, Mzo.*RPMs'*2*pi./60,'o')
xlabel('RPM')
ylabel('Power, W')