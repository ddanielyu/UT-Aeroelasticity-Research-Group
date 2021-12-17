function [Averages] = runSteady(StreamData,AvgData,steady_test)
%{

This function takes tripod 2021-2022 loadsprocessed data from runTripod.m 
and extracts the steady-only data into an Averages struct

Averages:
Cts_avg, Cts_err
Cps_avg, Cps_err
T_avg, T_err
P_avg, P_err
omega [rad/s], motor
sigma
rho
R
%}

%Run remaining LoadsProcessing.m Code (modified for Post-Aug 2021 Tripod Stand)
for i = 1:length(steady_test)
    %Find steady file in StreamData
    idx_steady = cellfun(@(x) strcmp(x, steady_test{i}), StreamData.names, 'UniformOutput', 1);

    %Create Averages structure
    Averages.cts_avg{i} = AvgData.avg_cts_inner{idx_steady};
    Averages.cts_err{i} = AvgData.err_cts_inner{idx_steady};
    Averages.cps_avg{i} = AvgData.avg_cps_inner{idx_steady};
    Averages.cps_err{i} = AvgData.err_cps_inner{idx_steady};
    Averages.rho{i} = StreamData.rho{idx_steady};
    Averages.OMEGA{i} = mean(StreamData.OMEGA{idx_steady}(2:length(StreamData.OMEGA{idx_steady})));
    Averages.R = StreamData.R;
    Averages.sigma = StreamData.sigma;
    Averages.T_avg{i} = Averages.cts_avg{i}*Averages.sigma*pi*Averages.R^2*....
        Averages.rho{i}*(Averages.OMEGA{i}*Averages.R)^2;
    Averages.T_err{i} = Averages.cts_err{i}*Averages.sigma*pi*Averages.R^2*....
        Averages.rho{i}*(Averages.OMEGA{i}*Averages.R)^2;
    Averages.P_avg{i} = Averages.cps_avg{i}*Averages.sigma*pi*Averages.R^2*....
        Averages.rho{i}*(Averages.OMEGA{i}*Averages.R)^3;
    Averages.P_err{i} = Averages.cps_err{i}*Averages.sigma*pi*Averages.R^2*....
        Averages.rho{i}*(Averages.OMEGA{i}*Averages.R)^3;
    
end
fprintf('\nSteady Processing done.\n');

end

