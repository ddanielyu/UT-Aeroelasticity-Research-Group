RPMs = MeanData.RPMs;
for i = 1:length(RPMs)
Fzo(i) = mean([StreamData.Fz_outer{i}]);
Fzi(i) = mean([StreamData.Fz_inner{i}]);
Mzo(i) = mean([StreamData.Mz_outer{i}]);
Mzi(i) = mean([StreamData.Mz_inner{i}]);
end
figure(4)
hold on
plot(RPMs, -Fzi*0.224809,'o')
xlabel('RPM')
ylabel('Thrust, lbf')
figure(5)
hold on
plot(RPMs, Mzi.*RPMs'*2*pi./60,'o')
xlabel('RPM')
ylabel('Power, W')