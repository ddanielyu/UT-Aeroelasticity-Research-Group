function [PhaseSync] = runPhaseSync(StreamData,phaseSync_test,offset,Torque_Trig,GR)
%{

This function processes data from phase-sync tests from Tripod Stand
2021-2022 Testing and outputs all phase-sync data into PhaseSync struc

PhaseSync:
time
rev
Thrust (sets separated by columns); T_avg; T_err;
Angle error "; Ang_err_avg; Ang_err_err;
Speed "; Speed_avg; Speed_err; (motor speed!!)
Torque "; Torque_avg; Torque_err; (motor torque!!)

%}

%% Constants
Fb = 8; %[N]
Mb = 1.5; %[N-m]

I75_b = 0.75; %CR Magnetics Current Sensor Bias (1% error) [A]
I150_b = 1.5;

Az_b = 360/1024/2; %Encoder resolution 1024 bits / rev

Fs = 10e3; %Sample rate [Hz]
pre_trig = 0.25;
post_trig = 1;


%% Find New Stream Data
for i = 1:length(phaseSync_test) 
    %Find phase-sync file in StreamData
    idx_PS = cellfun(@(x) strcmp(x, phaseSync_test{i}), StreamData.names, 'UniformOutput', 1);
    
    StreamData.Mz_inner{idx_PS} = abs(StreamData.Mz_inner{idx_PS});
    
    Time     = 0:1/Fs:(length(StreamData.encoder{idx_PS})-1)/Fs;
    
    %find angle trigger and calculate the rev-averaged angle error
    est_avg_Mz  = abs(mean(StreamData.Mz_inner{idx_PS}));
    
    if offset < 0
        start    = find((abs(StreamData.Mz_inner{idx_PS}) < est_avg_Mz - Torque_Trig),1);          %first index at which trig occurs
    else
        start    = find((abs(StreamData.Mz_inner{idx_PS}) > est_avg_Mz + Torque_Trig),1);          %first index at which trig occurs
    end

    rpm_motor0  = mean(StreamData.rpm{idx_PS}(2:start)*GR);                         %avg motor rpm prior to trig
    dazds       = rpm_motor0/60*360/Fs;                                     %rate of change in az over samples
    az_trig     = dazds*start;                                              %estimated az at the trig location [deg]
    ref_angle{idx_PS}   = dazds*(0:length(StreamData.encoder{idx_PS})-1);                                            %time array of entire data set up to trigger
    ref_angle{idx_PS}(start:end) = ref_angle{idx_PS}(start:end) + offset;

    %unwrap encoder measurements
    StreamData.unwrap_enc{idx_PS}(1) = 0;
    for ii = 2:length(StreamData.encoder{idx_PS}) 
        if StreamData.encoder{idx_PS}(ii) > StreamData.encoder{idx_PS}(ii-1)
            StreamData.unwrap_enc{idx_PS}(ii) = StreamData.unwrap_enc{idx_PS}(ii-1) + (StreamData.encoder{idx_PS}(ii) - StreamData.encoder{idx_PS}(ii-1));
        elseif StreamData.encoder{idx_PS}(ii) < StreamData.encoder{idx_PS}(ii-1)
            StreamData.unwrap_enc{idx_PS}(ii) = StreamData.unwrap_enc{idx_PS}(ii-1) + (StreamData.encoder{idx_PS}(ii) + (360 - StreamData.encoder{idx_PS}(ii-1)));
        end
    end

    %calc angle err
    idx = 1;
    while idx <= length(ref_angle{idx_PS})
        StreamData.angle_err{idx_PS}(idx) = StreamData.unwrap_enc{idx_PS}(idx) - ref_angle{idx_PS}(idx);
        idx = idx+1;
    end
    hi = StreamData.unwrap_enc{idx_PS} - StreamData.unwrap_enc{idx_PS}(start); 
    
    idx_trig        = (start - pre_trig*Fs):(start + post_trig*Fs); %index .25s pre-trig and 1s post-trig
    time            = Time - Time(start); time = time(idx_trig);%get time vect
    time0           = Time - Time(start); time0 = time0(idx_trig - 1);
    Mz_trans(:,idx_PS)   = StreamData.Mz_inner{idx_PS}(idx_trig);
    RPM_dazdt(:,idx_PS)  = (StreamData.unwrap_enc{idx_PS}(idx_trig) - StreamData.unwrap_enc{idx_PS}(idx_trig - 1))./(time - time0)*60/360;
    RPM_dazdt(:,idx_PS)  = fcleanup(RPM_dazdt(:,idx_PS), 'smoothdata', 'loess', 700);
    Angle_err(:,idx_PS)  = StreamData.angle_err{idx_PS}(idx_trig);
    Meas_angle(:,idx_PS) = hi(idx_trig);
    T_trans(:,idx_PS)   = abs(StreamData.Fz_inner{idx_PS}(idx_trig));
    T_trans(:,idx_PS)  = fcleanup(T_trans(:,idx_PS), 'smoothdata', 'loess', 1000);

    PhaseSync.time = time; %get time vect
    PhaseSync.Thrust(:,i) = T_trans(:,idx_PS);
    PhaseSync.Angle_err(:,i) = StreamData.angle_err{idx_PS}(idx_trig);
    PhaseSync.Speed(:,i) = RPM_dazdt(:,idx_PS);
    PhaseSync.Torque(:,i) = Mz_trans(:,idx_PS)/GR;
end

PhaseSync.rev = mean(Meas_angle');
PhaseSync.T_avg = mean(PhaseSync.Thrust');
PhaseSync.T_err = tinv(.975,size(PhaseSync.Thrust,2)) * sqrt(std(PhaseSync.Thrust').^2 + Fb^2)/sqrt(size(PhaseSync.Thrust,2));
PhaseSync.Ang_err_avg = mean(PhaseSync.Angle_err');
PhaseSync.Ang_err_err = tinv(.975,size(PhaseSync.Angle_err,2)) * sqrt(std(PhaseSync.Angle_err').^2 + Az_b^2)/sqrt(size(PhaseSync.Angle_err,2));
PhaseSync.Speed_avg = mean(PhaseSync.Speed');
PhaseSync.Speed_err = tinv(.975,size(PhaseSync.Speed,2)) * std(PhaseSync.Speed')/sqrt(size(PhaseSync.Speed,2));
PhaseSync.Torque_avg = mean(PhaseSync.Torque');
PhaseSync.Torque_err = tinv(.975,size(PhaseSync.Torque,2)) * sqrt(std(PhaseSync.Torque').^2 + Mb^2)/sqrt(size(PhaseSync.Torque,2));

fprintf('\nPhase-sync processing done.\n')

end
