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
Kt = 0.675; %[N-m/Arms] Estimated torque constant based on 12/16/2021 tripod stand testing of Q-I Curve Fit

Fb = 8; %[N]
Mb = 1.5; %[N-m]

I75_b = 0.75; %CR Magnetics Current Sensor Bias (1% error) [A]
I150_b = 1.5;
IQ_b = sqrt(2/3)*sqrt(3*I150_b^2); %error propagation of CR_150A current sensors

Az_b = 360/1024/2; %Encoder resolution 1024 bits / rev

Fs = 10e3; %Sample rate [Hz]
pre_trig = 0.25;
post_trig = 1;


%% Find New Stream Data
cnt = 1; %counter for splitting phase sync sets letters into cells
ij = 1; %counter for indexing each test within a set

for i = 1:length(phaseSync_test)
    if i > 1 
        current_test_name = split(phaseSync_test{i},{'_','-'});
        fprintf('\n%s%s\n','Processing phase-sync test ',current_test_name{3})        
        prev_test_name = split(phaseSync_test{i-1},{'_','-'});
        current_letter = current_test_name{3}; 
        previous_letter = prev_test_name{3};
        if current_letter ~= previous_letter; cnt = cnt+1; ij = 1; end %add another test set
    end
    
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
    
    Ref_angle = StreamData.unwrap_enc{idx_PS} - dazds*(0:length(StreamData.encoder{idx_PS})-1);
    
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
    Q_est_trans(:,idx_PS)   = StreamData.IQ{idx_PS}(idx_trig)*Kt/sqrt(2); %estimated motor torque based on IQ and Kt
    RPM_dazdt(:,idx_PS)  = (StreamData.unwrap_enc{idx_PS}(idx_trig) - StreamData.unwrap_enc{idx_PS}(idx_trig - 1))./(time - time0)*60/360;
    RPM_dazdt(:,idx_PS)  = fcleanup(RPM_dazdt(:,idx_PS), 'smoothdata', 'loess', 700);
    Angle_err(:,idx_PS)  = StreamData.angle_err{idx_PS}(idx_trig);
    Meas_angle(:,idx_PS) = hi(idx_trig);
    T_trans(:,idx_PS)   = abs(StreamData.Fz_inner{idx_PS}(idx_trig));
    T_trans(:,idx_PS)  = fcleanup(T_trans(:,idx_PS), 'smoothdata', 'loess', 1000);

    PhaseSync.Meas_angle{cnt}(:,ij) = Meas_angle(:,idx_PS); %get time vec
    PhaseSync.ref_angle{cnt}(:,ij) = Ref_angle(idx_trig);
    PhaseSync.time{cnt} = time; %get time vect
    PhaseSync.Thrust{cnt}(:,ij) = T_trans(:,idx_PS);
    PhaseSync.Angle_err{cnt}(:,ij) = StreamData.angle_err{idx_PS}(idx_trig);
    PhaseSync.Speed{cnt}(:,ij) = RPM_dazdt(:,idx_PS);
    PhaseSync.Torque{cnt}(:,ij) = Mz_trans(:,idx_PS)/GR;
    PhaseSync.Q_est{cnt}(:,ij) = Q_est_trans(:,idx_PS);
    [PhaseSync.Curr1_pk{cnt}(ij),curr_max_idx] = max(abs(StreamData.curr1{idx_PS}));
    PhaseSync.Torque_pk{cnt}(ij) = max(abs(StreamData.Mz_inner{idx_PS}))/GR;
    ij = ij + 1;
end

for i = 1:length(PhaseSync.Thrust)
    PhaseSync.rev{i} = mean(PhaseSync.Meas_angle{i}');
    PhaseSync.T_avg{i} = mean(PhaseSync.Thrust{i}');
    PhaseSync.T_err{i} = tinv(.975,size(PhaseSync.Thrust{i},2)) * sqrt(std(PhaseSync.Thrust{i}').^2 + Fb^2)/sqrt(size(PhaseSync.Thrust{i},2));
    PhaseSync.Ang_err_avg{i} = mean(PhaseSync.Angle_err{i}');
    PhaseSync.Ang_err_err{i} = tinv(.975,size(PhaseSync.Angle_err{i},2)) * sqrt(std(PhaseSync.Angle_err{i}').^2 + Az_b^2)/sqrt(size(PhaseSync.Angle_err{i},2));
    PhaseSync.ref_ang_avg{i} = mean(PhaseSync.ref_angle{i}');
    PhaseSync.ref_ang_err{i} = tinv(.975,size(PhaseSync.ref_angle{i},2)) * sqrt(std(PhaseSync.ref_angle{i}').^2 + Az_b^2)/sqrt(size(PhaseSync.ref_angle{i},2));
    PhaseSync.Speed_avg{i} = mean(PhaseSync.Speed{i}');
    PhaseSync.Speed_err{i} = tinv(.975,size(PhaseSync.Speed{i},2)) * std(PhaseSync.Speed{i}')/sqrt(size(PhaseSync.Speed{i},2));
    PhaseSync.Torque_avg{i} = mean(PhaseSync.Torque{i}');
    PhaseSync.Torque_err{i} = tinv(.975,size(PhaseSync.Torque{i},2)) * sqrt(std(PhaseSync.Torque{i}').^2 + Mb^2)/sqrt(size(PhaseSync.Torque{i},2));
    PhaseSync.Curr1_avg(i) = mean(PhaseSync.Curr1_pk{i});
    PhaseSync.Curr1_err(i) = tinv(.975,size(PhaseSync.Curr1_pk{i},2)) * sqrt(std(PhaseSync.Curr1_pk{i}).^2 + I150_b^2)/sqrt(size(PhaseSync.Curr1_pk{i},2));
    PhaseSync.Torque_pk_avg(i) = mean(PhaseSync.Torque_pk{i});
    PhaseSync.Torque_pk_err(i) = tinv(.975,size(PhaseSync.Torque_pk{i},2)) * sqrt(std(PhaseSync.Torque_pk{i}).^2 + Mb^2)/sqrt(size(PhaseSync.Torque_pk{i},2));
    PhaseSync.Q_est_avg{i} = mean(PhaseSync.Q_est{i}');
    PhaseSync.Q_est_err{i} = tinv(.975,size(PhaseSync.Q_est{i},2)) * sqrt(std(PhaseSync.Q_est{i}').^2 + (IQ_b*Kt/sqrt(2))^2)/sqrt(size(PhaseSync.Q_est{i},2));
end

fprintf('\nPhase-sync processing done.\n')

end

