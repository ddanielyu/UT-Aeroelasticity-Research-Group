function AvgData = fTotalAvg(RevData,SortedData,StreamData)
%{
EDITED ON: 03/15/2022
EDITED BY: MATT ASPER

Details: Dual Motor Processing Code. Accommodates 220315 Dual Motor Labview

%}
% CALCULATES TOTAL AVERAGE AND ERROR FOR STREAM DATA FILE
%
% INPUTS
%     RevData
%     StreamData
% OUTPUTS
%     AvgData      -> structure containting total average of nondimensional var in each streaming data file ( 1 value per file )
%                     .names            -> streaming data file names
%                     .avg_cts_outer    -> total average
%                     .avg_cps_outer
%                     .avg_cts_inner
%                     .avg_cps_outer
%                     .avg_FM_outer 
%                     .avg_FM_inner
%                     .avg_FM_tot 
%                     .avg_IQ 
%                     .err_cts_outer = sqrt(std(means)^2 + bias^2)/
%                                      sqrt(Nrev)
%                     .err_cps_outer   bias from uncertainty in load cell (8 N force, .48 N-m in moment; 
%                     .err_cts_inner
%                     .err_cps_outer
%                     .err_FM_outer 
%                     .err_FM_inner
%                     .err_FM_tot 
%                     .err_IQ
fprintf('\n%s\n', 'Averaging data');

for k = 1:length(StreamData.names)
    fprintf('\t%s', ['- ' SortedData.names{k} ' ... ']);
    
    OMEGA1 = nanmean(StreamData.OMEGA1{k});
    OMEGA2 = nanmean(StreamData.OMEGA2{k});

    % non-dimensionalization factor for CT
    ct_den1 = StreamData.rho{k} * (pi * StreamData.R^2) * (OMEGA1*StreamData.R).^2;
    ct_den2 = StreamData.rho{k} * (pi * StreamData.R^2) * (OMEGA2*StreamData.R).^2;

    % non-dimensionalization factor for CP
    cq_den1 = StreamData.rho{k} * (pi * StreamData.R^2) * (OMEGA1*StreamData.R).^2 * StreamData.R;
    cq_den2 = StreamData.rho{k} * (pi * StreamData.R^2) * (OMEGA2*StreamData.R).^2 * StreamData.R;

    AvgData.names{k} = StreamData.names{k};
    AvgData.avg_cts_outer{k} = nanmean(RevData.avg_cts_outer{k});
    AvgData.avg_cps_outer{k} = nanmean(RevData.avg_cps_outer{k});
    AvgData.avg_cts_inner{k} = nanmean(RevData.avg_cts_inner{k});
    AvgData.avg_cps_inner{k} = nanmean(RevData.avg_cps_inner{k});
    AvgData.avg_FM_outer{k} = nanmean(RevData.avg_FM_outer{k});
    AvgData.avg_FM_inner{k} = nanmean(RevData.avg_FM_inner{k});
    AvgData.avg_FM_tot{k} = nanmean(RevData.avg_FM_tot{k});
    AvgData.avg_ctcp{k} = nanmean(RevData.avg_ctcp{k});
    AvgData.avg_IQ1{k} = nanmean(RevData.avg_IQ1{k});
    AvgData.avg_IQ2{k} = nanmean(RevData.avg_IQ2{k});
    
    AvgData.avg_index{k} = nanmean(RevData.avg_index{k});

    cts_bias1 = 18.75 / ct_den1 / StreamData.sigma;
    cps_bias1 = 1.5 / cq_den1 / StreamData.sigma;
    cts_bias2 = 18.75 / ct_den2 / StreamData.sigma;
    cps_bias2 = 1.5 / cq_den2 / StreamData.sigma;
    
    I_bias = 1.5; %standard 1% error of 150 A sensor
    IQ_bias = sqrt(2/3)*sqrt(3*I_bias^2); %error propagation of CR_150A current sensors
    az_bias = 360/1024/2; %1024 bits/rev--> standard error half of smallest unit
    
    % standard error of the mean
    % se = std(means)/sqrt(Nrev)
    % multiply by 1.96 to get 95% CI
    AvgData.err_cts_outer{k} = 1.96* sqrt( std(RevData.ms_cts_outer{k})^2 + cts_bias2^2 )...
        / sqrt(SortedData.nrev2{k});
    AvgData.err_cps_outer{k} = 1.96* sqrt( std(RevData.ms_cps_outer{k})^2 + cps_bias2^2 )...
        / sqrt(SortedData.nrev2{k});
    AvgData.err_cts_inner{k} = 1.96* sqrt( std(RevData.ms_cts_inner{k})^2 + cts_bias1^2 )...
        / sqrt(SortedData.nrev1{k});
    AvgData.err_cps_inner{k} = 1.96* sqrt( std(RevData.ms_cps_inner{k})^2 + cps_bias1^2 )...
        / sqrt(SortedData.nrev1{k});

    AvgData.err_FM_outer{k} = 1.96* std(RevData.ms_FM_outer{k})/sqrt(SortedData.nrev2{k});
    AvgData.err_FM_inner{k} = 1.96* std(RevData.ms_FM_inner{k})/sqrt(SortedData.nrev1{k});
    AvgData.err_FM_tot{k} = 1.96* std(RevData.ms_FM_tot{k})/sqrt(SortedData.nrev1{k});
    
    AvgData.err_ctcp{k} = 1.96* std(RevData.ms_ctcp{k})/sqrt(SortedData.nrev1{k});
 
    AvgData.avg_cts_total{k} = (AvgData.avg_cts_inner{k} + AvgData.avg_cts_outer{k})/2;
    AvgData.avg_cps_total{k} = (AvgData.avg_cps_inner{k} + AvgData.avg_cps_outer{k})/2;
       
    AvgData.err_cts_total{k} = 1.96* sqrt( std(RevData.ms_cts_outer{k}+RevData.ms_cts_inner{k})^2 ...
    + (cts_bias1)^2 + (cts_bias2)^2 ) / sqrt(SortedData.nrev1{k}); 
    
    AvgData.err_cps_total{k} = 1.96* sqrt( std(RevData.ms_cps_outer{k}+RevData.ms_cps_inner{k})^2 ...
    + (cps_bias1)^2 + (cps_bias2)^2 ) / sqrt(SortedData.nrev1{k});

    AvgData.err_IQ1{k} = 1.96* sqrt( std(RevData.ms_IQ1{k})^2 + IQ_bias^2)/sqrt(SortedData.nrev1{k});
    AvgData.err_IQ2{k} = 1.96* sqrt( std(RevData.ms_IQ2{k})^2 + IQ_bias^2)/sqrt(SortedData.nrev2{k});
    
    AvgData.err_index{k} = 1.96* sqrt( std(RevData.ms_index{k})^2 + az_bias^2)/sqrt(SortedData.nrev1{k});
    
    fprintf('%s\n', ' Ok');
end


