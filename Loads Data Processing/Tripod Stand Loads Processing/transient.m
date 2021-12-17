function [] = transient
%{
This function processes phase-sync transient data from the Fall 2021 Tripod
Stand Testing

Written By: Matt Asper
Date: 01 Dec 2021

%}

%% Constants
Fb = 18.75; %[N]
Mb = 1.5; %[N-m]
Fs = 10e3;
GR = 1.2;

I75_b = 0.75; %CR Magnetics Current Sensor Bias (1% error) [A]
I150_b = 1.5;

Az_b = 360/1024/2; %Encoder resolution 1024 bits / rev

Fs = 10e3; %Sample rate [Hz]
pre_trig = 0.25;
post_trig = 1;
offset = 5; %angle offset [deg]

%% Load Data
fprintf('Select Data File to Load...\n');
[data_file,data_path] = uigetfile();
load(fullfile(data_path,data_file),'StreamData');

Pos_Angles_idx = input('Positive Offset [indicies]: ');
Neg_Angles_idx = input('Negative Offset [indicies]: ');


if isempty(Pos_Angles_idx) ~= 1
    Pos_Mz = StreamData.Mz_inner{Pos_Angles_idx};
    pos_enc = StreamData.encoder{Pos_Angles_idx};
    pos_rpm = StreamData.rpm{Pos_Angles_idx};
end

if isempty(Pos_Angles_idx) ~= 1
    Neg_Mz = StreamData.Mz_inner{Neg_Angles_idx};
    neg_enc = StreamData.encoder{Neg_Angles_idx};
    neg_rpm = StreamData.rpm{Neg_Angles_idx};
end

%% Process Data

for i = 1:length(Pos_Mz)
    Pos_Mz{i} = abs(Pos_Mz{i});
    
    Time     = 0:1/Fs:(length(Pos_Mz{i})-1)/Fs;

    %find angle trigger and calculate the rev-averaged angle error
    est_avg_Mz  = mean(Pos_Mz{i});

    start    = find((Pos_Mz{i} > est_avg_Mz + 5),1);          %first index at which trig occurs

    rpm_motor0  = mean(pos_rpm{i}*GR);                         %avg motor rpm prior to trig
    dazds       = rpm_motor0/60*360/Fs;                                     %rate of change in az over samples
    az_trig     = dazds*start;                                              %estimated az at the trig location [deg]
    ref_angle{i}   = dazds*(0:length(pos_enc{i})-1);                                            %time array of entire data set up to trigger
    ref_angle{i}(start:end) = ref_angle{i}(start:end) + offset;

    %unwrap encoder measurements
    StreamData.unwrap_enc{i}(1) = 0;
    for ii = 2:length(StreamData.encoder{i}) 
        if StreamData.encoder{i}(ii) > StreamData.encoder{i}(ii-1)
            StreamData.unwrap_enc{i}(ii) = StreamData.unwrap_enc{i}(ii-1) + (StreamData.encoder{i}(ii) - StreamData.encoder{i}(ii-1));
        elseif StreamData.encoder{i}(ii) < StreamData.encoder{i}(ii-1)
            StreamData.unwrap_enc{i}(ii) = StreamData.unwrap_enc{i}(ii-1) + (StreamData.encoder{i}(ii) + (360 - StreamData.encoder{i}(ii-1)));
        end
    end

    %calc angle err
    idx = 1;
    while idx <= length(ref_angle{i})
        StreamData.angle_err{i}(idx) = StreamData.unwrap_enc{i}(idx) - ref_angle{i}(idx);
        idx = idx+1;
    end
    hi = StreamData.unwrap_enc{i} - StreamData.unwrap_enc{i}(start(i)); 

    idx_trig        = (start(i) - pre_trig*Fs):(start(i) + post_trig*Fs); %index .25s pre-trig and 1s post-trig
    time            = Time - Time(start(i)); time = time(idx_trig);%get time vect
    time0           = Time - Time(start(i)); time0 = time0(idx_trig - 1);
    Mz_trans(:,i)   = StreamData.Mz_inner{i}(idx_trig);
    RPM_dazdt(:,i)  = (StreamData.unwrap_enc{i}(idx_trig) - StreamData.unwrap_enc{i}(idx_trig - 1))./(time - time0)*60/360;
    RPM_dazdt(:,i)  = fcleanup(RPM_dazdt(:,i), 'smoothdata', 'loess', 700);
    Angle_err(:,i)  = StreamData.angle_err{i}(idx_trig);
    Meas_angle(:,i) = hi(idx_trig);
    
end

avg_Mz_trans = mean(Mz_trans');
err_Mz_trans = 2.228 * sqrt(std(Mz_trans').^2 + Mb^2)/sqrt(10);
avg_angle_err = mean(angle_err');
err_angle_err = 2.228 * sqrt(std(angle_err').^2 + Az_b^2)/sqrt(10);
Meas_angle    = mean(Meas_angle');
avg_RPM_dazdt = mean(RPM_dazdt');
err_RPM_dazdt = 2.228 * std(RPM_dazdt')/sqrt(10);

end