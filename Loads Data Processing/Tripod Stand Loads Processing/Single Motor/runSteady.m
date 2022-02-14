function [Averages] = runSteady(StreamData,AvgData,steady_test)
%{
EDITED ON: 01/12/2022
EDITED BY: MATT ASPER

Details: Includes estimated torque based on IQ and Kt of the EMRAX 208

This function takes tripod 2021-2022 loadsprocessed data from runTripod.m 
and extracts the steady-only data into an Averages struct

Averages:
Cts_avg, Cts_err
Cps_avg, Cps_err
T_avg, T_err
P_avg, P_err
Q_est_avg, Q_est_err [estimated MOTOR torque]
omega [rad/s], motor
sigma
rho
R
%}

Kt = 0.675; %Estimated torque constant based on 12/16/2021 tripod stand testing of Q-I Curve Fit


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
    Averages.Q_est_avg{i} = AvgData.avg_IQ{idx_steady}*Kt;
    Averages.Q_est_err{i} = AvgData.err_IQ{idx_steady}*Kt;
    
end
fprintf('\nSteady Processing done.\n');

end

