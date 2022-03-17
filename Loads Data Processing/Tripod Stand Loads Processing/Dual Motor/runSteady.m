function [Averages] = runSteady(StreamData,AvgData,steady_test)
%{
EDITED ON: 03/15/2022
EDITED BY: MATT ASPER

Details: Dual Motor Processing Code. Accommodates 220315 Dual Motor Labview

Averages:
Cts_avg, Cts_err
Cps_avg, Cps_err
T_avg, T_err
P_avg, P_err
FM_up,FM_lo,FM_tot
P1_est_avg, P1_est_err [estimated servo ROTOR Power]
P2_est_avg, P2_est_err [estimated follower ROTOR Power]
omega1 [rad/s], servo rotor
omega2 [rad/s], follower rotor
index
sigma
rho
R
%}

Kt = 0.675; %Estimated torque constant based on 12/16/2021 tripod stand testing of Q-I Curve Fit


%Run remaining LoadsProcessing.m Code (modified for Post-Aug 2021 Tripod Stand)
for i = 1:length(steady_test)
    %Find steady file in StreamData
    idx_steady = cellfun(@(x) strcmp(x, steady_test{i}), StreamData.names, 'UniformOutput', 1);

    Averages.rho{i} = StreamData.rho{idx_steady};
    Averages.OMEGA1{i} = mean(StreamData.OMEGA1{idx_steady}(2:length(StreamData.OMEGA1{idx_steady})));
    Averages.OMEGA2{i} = mean(StreamData.OMEGA2{idx_steady}(2:length(StreamData.OMEGA2{idx_steady})));
    Averages.R = StreamData.R;
    Averages.sigma = StreamData.sigma;
    Averages.GR = StreamData.GR;
    
    %Create Averages structure
    %upper
    Averages.cts_up_avg{i} = AvgData.avg_cts_inner{idx_steady};
    Averages.cts_up_err{i} = AvgData.err_cts_inner{idx_steady};
    Averages.cps_up_avg{i} = AvgData.avg_cps_inner{idx_steady};
    Averages.cps_up_err{i} = AvgData.err_cps_inner{idx_steady};
    Averages.cps_up_est_avg{i} = AvgData.avg_IQ1{idx_steady}*Kt/sqrt(2)*Averages.OMEGA1{i}*Averages.GR/...
        (Averages.sigma*pi*Averages.R^2*Averages.rho{i}*(Averages.OMEGA1{i}*Averages.R)^3);
    Averages.cps_up_est_err{i} = AvgData.err_IQ1{idx_steady}*Kt/sqrt(2)*Averages.OMEGA1{i}*Averages.GR/...
        (Averages.sigma*pi*Averages.R^2*Averages.rho{i}*(Averages.OMEGA1{i}*Averages.R)^3);
    
    %lower
    Averages.cts_lo_avg{i} = AvgData.avg_cts_outer{idx_steady};
    Averages.cts_lo_err{i} = AvgData.err_cts_outer{idx_steady};
    Averages.cps_lo_avg{i} = AvgData.avg_cps_outer{idx_steady};
    Averages.cps_lo_err{i} = AvgData.err_cps_outer{idx_steady};
    Averages.cps_lo_est_avg{i} = AvgData.avg_IQ2{idx_steady}*Kt/sqrt(2)*Averages.OMEGA2{i}*Averages.GR/...
        (Averages.sigma*pi*Averages.R^2*Averages.rho{i}*(Averages.OMEGA2{i}*Averages.R)^3);
    Averages.cps_lo_est_err{i} = AvgData.err_IQ2{idx_steady}*Kt/sqrt(2)*Averages.OMEGA2{i}*Averages.GR/...
        (Averages.sigma*pi*Averages.R^2*Averages.rho{i}*(Averages.OMEGA2{i}*Averages.R)^3);
    
    %total
    Averages.cts_avg{i} = mean([Averages.cts_up_avg{i},Averages.cts_lo_avg{i}]);
    Averages.cts_err{i} = mean([Averages.cts_up_err{i},Averages.cts_lo_err{i}]);
    Averages.cps_avg{i} = mean([Averages.cps_up_avg{i},Averages.cps_lo_avg{i}]);
    Averages.cps_err{i} = mean([Averages.cps_up_err{i},Averages.cps_lo_err{i}]);
    Averages.cps_est_avg{i} = mean([Averages.cps_up_est_avg{i},Averages.cps_lo_est_avg{i}]);
    Averages.cps_est_err{i} = mean([Averages.cps_up_est_err{i},Averages.cps_lo_est_err{i}]);
    
    %FM
    Averages.FM_lo{i} = AvgData.avg_FM_outer{idx_steady}; Averages.FM_err_lo{i} = AvgData.err_FM_outer{idx_steady};
    Averages.FM_up{i} = AvgData.avg_FM_inner{idx_steady}; Averages.FM_err_up{i} = AvgData.err_FM_inner{idx_steady};
    Averages.FM_tot{i}= AvgData.avg_FM_tot{idx_steady}; Averages.FM_err_tot{i} = AvgData.err_FM_tot{idx_steady};
    
    %upper
    Averages.T_avg_up{i} = Averages.cts_up_avg{i}*Averages.sigma*pi*Averages.R^2*....
        Averages.rho{i}*(Averages.OMEGA1{i}*Averages.R)^2;
    Averages.T_err_up{i} = Averages.cts_up_err{i}*Averages.sigma*pi*Averages.R^2*....
        Averages.rho{i}*(Averages.OMEGA1{i}*Averages.R)^2;
    Averages.P_avg_up{i} = Averages.cps_up_avg{i}*Averages.sigma*pi*Averages.R^2*....
        Averages.rho{i}*(Averages.OMEGA1{i}*Averages.R)^3;
    Averages.P_err_up{i} = Averages.cps_up_err{i}*Averages.sigma*pi*Averages.R^2*....
        Averages.rho{i}*(Averages.OMEGA1{i}*Averages.R)^3;
    Averages.P_est_avg_up{i} = AvgData.avg_IQ1{idx_steady}*Kt/sqrt(2)*Averages.OMEGA1{i}*Averages.GR;
    Averages.P_est_err_up{i} = AvgData.err_IQ1{idx_steady}*Kt/sqrt(2)*Averages.OMEGA1{i}*Averages.GR;
    
    %lower
    Averages.T_avg_lo{i} = Averages.cts_lo_avg{i}*Averages.sigma*pi*Averages.R^2*....
        Averages.rho{i}*(Averages.OMEGA2{i}*Averages.R)^2;
    Averages.T_err_lo{i} = Averages.cts_lo_err{i}*Averages.sigma*pi*Averages.R^2*....
        Averages.rho{i}*(Averages.OMEGA2{i}*Averages.R)^2;
    Averages.P_avg_lo{i} = Averages.cps_lo_avg{i}*Averages.sigma*pi*Averages.R^2*....
        Averages.rho{i}*(Averages.OMEGA2{i}*Averages.R)^3;
    Averages.P_err_lo{i} = Averages.cps_lo_err{i}*Averages.sigma*pi*Averages.R^2*....
        Averages.rho{i}*(Averages.OMEGA2{i}*Averages.R)^3;
    Averages.P_est_avg_lo{i} = AvgData.avg_IQ2{idx_steady}*Kt/sqrt(2)*Averages.OMEGA2{i}*Averages.GR;
    Averages.P_est_err_lo{i} = AvgData.err_IQ2{idx_steady}*Kt/sqrt(2)*Averages.OMEGA2{i}*Averages.GR;
    
    
    %total
    Averages.T_avg{i} = Averages.T_avg_lo{i} + Averages.T_avg_up{i};
    Averages.T_err{i} = Averages.T_err_lo{i} + Averages.T_err_up{i};
    Averages.P_avg{i} = Averages.P_avg_lo{i} + Averages.P_avg_up{i};
    Averages.P_err{i} = Averages.P_err_lo{i} + Averages.P_err_up{i};
    Averages.P_est_avg{i} = Averages.P_est_avg_lo{i} + Averages.P_est_avg_up{i};
    Averages.P_est_err{i} = Averages.P_est_err_lo{i} + Averages.P_est_err_up{i};
    
    %index
    Averages.index_avg{i} = AvgData.avg_index{idx_steady};
    Averages.index_err{i} = AvgData.err_index{idx_steady};
end
fprintf('\nSteady Processing done.\n');

end

