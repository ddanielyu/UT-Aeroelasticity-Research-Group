function [StreamData, SortedData] = fSortStream(StreamData)
%{
EDITED ON: 03/15/2022
EDITED BY: MATT ASPER

Details: Goes with Dual Motor Tripod Labview_rev 220315


%}
%INPUTS
%     StreamData
%     conditions = [Temperature [F], % Humidity, Prassure [in-hg]]
% OUTPUTS
%     StreamData -> adds calculated variables to StreamData
%                     .sigma
%                     .rho
%                     .OMEGA
%                     .binsize
%     SortedData -> structure containing streaming data sorted into matrices;
%                   each cell = one mean data file; each row matrix = one revolution
%                   calculates ct/sigma,cp/sigma, and FM's
%                     .binsize      -> number of points recorded in each revolution
%                     .nrevs1
%                     .nrevs2
%                     .names
%                     .check
%                     .encoder1      -> servo azimuthal position
%                     .azimuth1      -> vector of azimuths at which all data
%                                      is resampled
%                     .encoder2      -> follower azimuthal position
%                     .azimuth2      -> vector of azimuths at which all data
%                                      is resampled
%                     .instRPM1      -> instantaneous servo rotor speed at each time
%                     .instRPM2      -> instantaneous follower rotor speed at each time
%                     .Fx_outer     -> all data points recorded by labview
%                     .Fy_outer
%                     .Fz_outer
%                     .Mx_outer
%                     .My_outer
%                     .Mz_outer
%                     .Fx_inner
%                     .Fy_inner
%                     .Fz_inner
%                     .Mx_inner
%                     .My_inner
%                     .Mz_inner
%                     .cts_outer
%                     .cps_outer
%                     .cts_inner
%                     .cps_outer
%                     .FM_outer = (F)^1.5 / P / sqrt(2 rho A)
%                     .FM_inner
%                     .FM_tot = (F + F)^1.5 / (P+P) / sqrt(2 rho A)

%% DEFINE CONSTANTS
Nbits = 12;       % number of bits in the 
StreamData.R = 1.016;
c = 0.08;
Nblade = 2;
StreamData.sigma = StreamData.R*c*Nblade / (pi*StreamData.R^2);

SR = 4000; % SAMPLE RATE
Naz = 1000;   % dpsi = 0.36 deg
enc = 'y';

%% CALCULATE OMEGA FROM 1/Rev
for k = 1:length(StreamData.names)
    
    StreamData.binsize1{k} = zeros(1,StreamData.nrev1{k});
    StreamData.binsize2{k} = zeros(1,StreamData.nrev2{k});
    
    for n = 1:StreamData.nrev1{k}
        StreamData.binsize1{k}(n) = sum(StreamData.revolution1{k}(:) == n-1);
    end

    for n = 1:StreamData.nrev2{k}
        StreamData.binsize2{k}(n) = sum(StreamData.revolution2{k}(:) == n-1);
    end
            
    
    StreamData.OMEGA1{k} = SR./StreamData.binsize1{k} * 2 * pi;
    SortedData.binsize1{k} = StreamData.binsize1{k};
    SortedData.nrev1{k} = StreamData.nrev1{k};
    
    StreamData.OMEGA2{k} = SR./StreamData.binsize2{k} * 2 * pi;
    SortedData.binsize2{k} = StreamData.binsize2{k};
    SortedData.nrev2{k} = StreamData.nrev2{k};
    
    SortedData.names{k} = StreamData.names{k};
end

fprintf('\n%s\n', 'Sorting data');

%% CALCULATE CT/S AND CP/S
for k = 1:length(StreamData.names)
    
    fprintf('\t%s', ['- ' StreamData.names{k} ' ... ']);
        
    SortedData.check1{k} = [];
    SortedData.check2{k} = [];
    
    SortedData.encoder1{k} = [];
    SortedData.azimuth1{k} = [];
    
    SortedData.encoder2{k} = [];
    SortedData.azimuth2{k} = [];
    
    SortedData.index{k} = [];

    SortedData.instRPM1{k} = [];    
    SortedData.instRPM2{k} = [];    
    SortedData.Fx_outer{k} = [];
    SortedData.Fy_outer{k} = [];
    SortedData.Fz_outer{k} = []; 
    SortedData.Mx_outer{k} = [];
    SortedData.My_outer{k} = [];
    SortedData.Mz_outer{k} = [];
    SortedData.Fx_inner{k} = [];
    SortedData.Fy_inner{k} = [];
    SortedData.Fz_inner{k} = [];
    SortedData.Mx_inner{k} = [];
    SortedData.My_inner{k} = [];
    SortedData.Mz_inner{k} = [];

    SortedData.ax{k} = [];
    SortedData.ay{k} = [];
    
    SortedData.curr1{k} = [];
    SortedData.curr2{k} = [];
    SortedData.curr3{k} = [];
    SortedData.curr4{k} = [];
    SortedData.curr5{k} = [];
    SortedData.curr6{k} = [];
    SortedData.bus_curr{k} = [];
    SortedData.IQ1{k} = [];
    SortedData.IQ2{k} = [];
    
    SortedData.index{k} = [];
    
    SortedData.cts_outer{k} = [];
    SortedData.cps_outer{k} = [];
    SortedData.cts_inner{k} = [];
    SortedData.cps_inner{k} = [];
    
    SortedData.FM_outer{k} = [];
    SortedData.FM_inner{k} = [];
    SortedData.FM_tot{k} = [];
    
    SortedData.ctcp{k} = [];
    
    
    SortedData.azimuth1{k} = linspace(0, 360*(1-1/Naz), Naz);
    SortedData.azimuth2{k} = linspace(0, 360*(1-1/Naz), Naz);
    
    %% Servo Loop (Inner)
    count1 = SortedData.binsize1{k}(1)+1;
    
    for n = 2:StreamData.nrev1{k} %cut out first rev
        b = SortedData.binsize1{k}(n);
        
        SortedData.check1{k}(n,1:b) = StreamData.revolution1{k}(count1:count1-1+b)';
        
        %"zero" enc for each rev.
        az = StreamData.encoder1{k}(count1:count1-1+b)' - StreamData.encoder1{k}(count1); 
        
        % SortedData.encoder{k}(n,1:b) = StreamData.encoder{k}(count:count-1+b)';
        % if binsize is greater than maximum readable by encoder, decimate
        % the data stream by dcm8 = 3
        dcm8 = 3;
        if b > (2^Nbits-1)
            az = az(1:dcm8:end);     % length of decimated az is ceil(b/3);
            SortedData.encoder1{k}(n,1:length(az)) = az;
            azdt = circshift(az,-1);
            instRPM = (azdt(1:end-1) - az(1:end-1)) *SR * pi /180 /dcm8; % instantaneous RPM, rad/s
            instRPM = [instRPM instRPM(end)]; % add one element to get size 1xb
            % interpolate to azimuth with dpsi = 1/Naz
            SortedData.instRPM1{k}(n,:) = interp1(az, instRPM, SortedData.azimuth1{k}, 'pchip');

            SortedData.Fx_inner{k}(n,:) = interp1(az, StreamData.Fx_inner{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.Fy_inner{k}(n,:) = interp1(az, StreamData.Fy_inner{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.Fz_inner{k}(n,:) = interp1(az, StreamData.Fz_inner{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.Mx_inner{k}(n,:) = interp1(az, StreamData.Mx_inner{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.My_inner{k}(n,:) = interp1(az, StreamData.My_inner{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.Mz_inner{k}(n,:) = interp1(az, StreamData.Mz_inner{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.ax{k}(n,:) = interp1(az, StreamData.ax{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.ay{k}(n,:) = interp1(az, StreamData.ay{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.curr1{k}(n,:) = interp1(az, StreamData.curr1{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.curr2{k}(n,:) = interp1(az, StreamData.curr2{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.curr3{k}(n,:) = interp1(az, StreamData.curr3{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.bus_curr{k}(n,:) = interp1(az, StreamData.bus_curr{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.IQ1{k}(n,:) = interp1(az, StreamData.IQ1{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
            
            SortedData.index{k}(n,:) = interp1(az, StreamData.index{k}(count1:dcm8:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
            
        else
            SortedData.encoder1{k}(n,1:b) = az;
            azdt = circshift(az,-1);
            instRPM = (azdt(1:end-1) - az(1:end-1)) *SR * pi /180; % instantaneous RPM, rad/s
            instRPM = [instRPM instRPM(end)]; % add one element to get size 1xb
            % interpolate to azimuth with dpsi = 1/Naz

            SortedData.instRPM1{k}(n,:) = interp1(az, instRPM, SortedData.azimuth1{k}, 'pchip');

            SortedData.Fx_inner{k}(n,:) = interp1(az, StreamData.Fx_inner{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.Fy_inner{k}(n,:) = interp1(az, StreamData.Fy_inner{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.Fz_inner{k}(n,:) = interp1(az, StreamData.Fz_inner{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.Mx_inner{k}(n,:) = interp1(az, StreamData.Mx_inner{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.My_inner{k}(n,:) = interp1(az, StreamData.My_inner{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.Mz_inner{k}(n,:) = interp1(az, StreamData.Mz_inner{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.ax{k}(n,:) = interp1(az, StreamData.ax{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.ay{k}(n,:) = interp1(az, StreamData.ay{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.curr1{k}(n,:) = interp1(az, StreamData.curr1{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.curr2{k}(n,:) = interp1(az, StreamData.curr2{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.curr3{k}(n,:) = interp1(az, StreamData.curr3{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.bus_curr{k}(n,:) = interp1(az, StreamData.bus_curr{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

            SortedData.IQ1{k}(n,:) = interp1(az, StreamData.IQ1{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
            
            SortedData.index{k}(n,:) = interp1(az, StreamData.index{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
        end
        
        count1 = count1+b;
                
    end
    
    
    %account for mis-match in revs
%     if SortedData.nrev1{k} > SortedData.nrev2{k}
%         SortedData.instRPM1{k}(end,:) = [];
% 
%         SortedData.Fx_inner{k}(end,:) = interp1(az, StreamData.Fx_inner{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.Fy_inner{k}(end,:) = interp1(az, StreamData.Fy_inner{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.Fz_inner{k}(end,:) = interp1(az, StreamData.Fz_inner{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.Mx_inner{k}(end,:) = interp1(az, StreamData.Mx_inner{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.My_inner{k}(end,:) = interp1(az, StreamData.My_inner{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.Mz_inner{k}(end,:) = interp1(az, StreamData.Mz_inner{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.ax{k}(end,:) = interp1(az, StreamData.ax{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.ay{k}(end,:) = interp1(az, StreamData.ay{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.curr1{k}(,:) = interp1(az, StreamData.curr1{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.curr2{k}(n,:) = interp1(az, StreamData.curr2{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.curr3{k}(n,:) = interp1(az, StreamData.curr3{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.bus_curr{k}(n,:) = interp1(az, StreamData.bus_curr{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.IQ1{k}(n,:) = interp1(az, StreamData.IQ1{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');
% 
%         SortedData.index{k}(n,:) = interp1(az, StreamData.index{k}(count1:count1-1+b)', SortedData.azimuth1{k}, 'pchip');

        

%% Follower Loop

    count2 = SortedData.binsize2{k}(1)+1;
    
    for n = 2:StreamData.nrev2{k} %cut out first rev
        b = SortedData.binsize2{k}(n);
        
        SortedData.check2{k}(n,1:b) = StreamData.revolution2{k}(count2:count2-1+b)';

        %"zero" enc for each rev.
        az = StreamData.encoder2{k}(count2:count2-1+b)' - StreamData.encoder2{k}(count2);

        % SortedData.encoder{k}(n,1:b) = StreamData.encoder{k}(count:count-1+b)';
        % if binsize is greater than maximum readable by encoder, decimate
        % the data stream by dcm8 = 3
        dcm8 = 3;
        if b > (2^Nbits-1)
            az = az(1:dcm8:end);     % length of decimated az is ceil(b/3);
            SortedData.encoder2{k}(n,1:length(az)) = az;
            azdt = circshift(az,-1);
            instRPM = (azdt(1:end-1) - az(1:end-1)) *SR * pi /180 /dcm8; % instantaneous RPM, rad/s
            instRPM = [instRPM instRPM(end)]; % add one element to get size 1xb
            
            % interpolate to azimuth with dpsi = 1/Naz
            SortedData.instRPM2{k}(n,:) = interp1(az, instRPM, SortedData.azimuth2{k}, 'pchip');

            SortedData.Fx_outer{k}(n,:) = interp1(az, StreamData.Fx_outer{k}(count2:dcm8:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.Fy_outer{k}(n,:) = interp1(az, StreamData.Fy_outer{k}(count2:dcm8:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.Fz_outer{k}(n,:) = interp1(az, StreamData.Fz_outer{k}(count2:dcm8:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.Mx_outer{k}(n,:) = interp1(az, StreamData.Mx_outer{k}(count2:dcm8:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.My_outer{k}(n,:) = interp1(az, StreamData.My_outer{k}(count2:dcm8:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.Mz_outer{k}(n,:) = interp1(az, StreamData.Mz_outer{k}(count2:dcm8:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.curr4{k}(n,:) = interp1(az, StreamData.curr4{k}(count2:dcm8:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.curr5{k}(n,:) = interp1(az, StreamData.curr5{k}(count2:dcm8:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.curr6{k}(n,:) = interp1(az, StreamData.curr5{k}(count2:dcm8:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.IQ2{k}(n,:) = interp1(az, StreamData.IQ2{k}(count2:dcm8:count2-1+b)', SortedData.azimuth2{k}, 'pchip');
        
        else
            SortedData.encoder2{k}(n,1:b) = az;
            azdt = circshift(az,-1);
            instRPM = (azdt(1:end-1) - az(1:end-1)) *SR * pi /180; % instantaneous RPM, rad/s
            instRPM = [instRPM instRPM(end)]; % add one element to get size 1xb
            % interpolate to azimuth with dpsi = 1/Naz

            SortedData.instRPM2{k}(n,:) = interp1(az, instRPM, SortedData.azimuth2{k}, 'pchip');


            SortedData.Fx_outer{k}(n,:) = interp1(az, StreamData.Fx_outer{k}(count2:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.Fy_outer{k}(n,:) = interp1(az, StreamData.Fy_outer{k}(count2:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.Fz_outer{k}(n,:) = interp1(az, StreamData.Fz_outer{k}(count2:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.Mx_outer{k}(n,:) = interp1(az, StreamData.Mx_outer{k}(count2:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.My_outer{k}(n,:) = interp1(az, StreamData.My_outer{k}(count2:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.Mz_outer{k}(n,:) = interp1(az, StreamData.Mz_outer{k}(count2:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.curr4{k}(n,:) = interp1(az, StreamData.curr4{k}(count2:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.curr5{k}(n,:) = interp1(az, StreamData.curr5{k}(count2:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.curr6{k}(n,:) = interp1(az, StreamData.curr6{k}(count2:count2-1+b)', SortedData.azimuth2{k}, 'pchip');

            SortedData.IQ2{k}(n,:) = interp1(az, StreamData.IQ2{k}(count2:count2-1+b)', SortedData.azimuth1{k}, 'pchip');
        end
        
        count2 = count2+b;
    end
    
            SortedData.Fx_outer{k}(1,:) = [];
            SortedData.Fy_outer{k}(1,:) = [];              
            SortedData.Fz_outer{k}(1,:) = [];                
            SortedData.Mx_outer{k}(1,:) = [];
            SortedData.My_outer{k}(1,:) = [];
            SortedData.Mz_outer{k}(1,:) = [];
            SortedData.Fx_inner{k}(1,:) = [];
            SortedData.Fy_inner{k}(1,:) = [];
            SortedData.Fz_inner{k}(1,:) = [];
            SortedData.Mx_inner{k}(1,:) = [];
            SortedData.My_inner{k}(1,:) = [];
            SortedData.Mz_inner{k}(1,:) = [];
            SortedData.ax{k}(1,:) = [];
            SortedData.ay{k}(1,:) = [];
            SortedData.curr1{k}(1,:) = [];
            SortedData.curr2{k}(1,:) = [];
            SortedData.curr3{k}(1,:) = [];
            SortedData.curr4{k}(1,:) = [];
            SortedData.curr5{k}(1,:) = [];
            SortedData.curr6{k}(1,:) = [];
            SortedData.bus_curr{k}(1,:) = [];
            SortedData.IQ1{k}(1,:) = [];
            SortedData.IQ2{k}(1,:) = [];
            SortedData.index{k}(1,:) = [];
    
    OMEGA1 = StreamData.OMEGA1{k}(2:length(StreamData.OMEGA1{k}));
    OMEGA2 = StreamData.OMEGA2{k}(2:length(StreamData.OMEGA2{k}));

    SortedData.cts_outer{k} = SortedData.Fz_outer{k} ./ StreamData.rho{k} / (pi * StreamData.R^2) ./ (OMEGA2'*StreamData.R).^2 / StreamData.sigma;
    SortedData.cps_outer{k} = SortedData.Mz_outer{k} ./ StreamData.rho{k} / (pi * StreamData.R^2) ./ (OMEGA2'*StreamData.R).^2 /StreamData.R / StreamData.sigma;
    SortedData.cts_inner{k} = SortedData.Fz_inner{k} ./ StreamData.rho{k} / (pi * StreamData.R^2) ./ (OMEGA1'*StreamData.R).^2 / StreamData.sigma;
    SortedData.cps_inner{k} = SortedData.Mz_inner{k} ./ StreamData.rho{k} / (pi * StreamData.R^2) ./ (OMEGA1'*StreamData.R).^2 /StreamData.R / StreamData.sigma;

    SortedData.FM_outer{k} = sqrt(StreamData.sigma) * (abs(SortedData.cts_outer{k}).^(3/2)) ./ sqrt(2)./ (SortedData.cps_outer{k});
    SortedData.FM_inner{k} = sqrt(StreamData.sigma) * (abs(SortedData.cts_inner{k}).^(3/2)) ./ sqrt(2)./ (SortedData.cps_inner{k});
    SortedData.FM_tot{k} = sqrt(StreamData.sigma) * ...
        (abs(SortedData.cts_outer{k} + SortedData.cts_inner{k}).^(3/2)) ./ sqrt(2)./ (SortedData.cps_outer{k} + SortedData.cps_inner{k});
    
    SortedData.ctcp{k} = (SortedData.cts_outer{k}+SortedData.cts_inner{k})./(SortedData.cps_outer{k}+SortedData.cps_inner{k});
    fprintf('%s\n', ['Servo/Follower Nrevs: ' num2str(StreamData.nrev1{k}) '/' num2str(StreamData.nrev2{k}) '... Ok'] );
   
end

end


