function [PhaseSync] = runPhaseSync(StreamData,phaseSync_test,offset,Trig,GR,slew,SR)
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
Kt = StreamData.Kt; %[N-m/Arms] Estimated torque constant based on 12/16/2021 tripod stand testing of Q-I Curve Fit

Fb = 18.75; %[N]
Mb = 1.5; %[N-m]

I75_b = 0.75; %CR Magnetics Current Sensor Bias (1% error) [A]
I150_b = 1.5;
IQ_b = sqrt(2/3)*I150_b; %error propagation of CR_150A current sensors

Az_b = 360/1024/2; %Encoder resolution 1024 bits / rev

pre_trig = 1;
post_trig = 4;

%% Find New Stream Data

for i = 1:length(phaseSync_test)

    Time     = 0:1/SR:(length(StreamData.encoder{i})-1)/SR;
    
    %estimate the avg rpm based on first .5s of data
    est_avg_rpm = abs(nanmean(StreamData.rpm{i}(1:SR/2)));
    
    if offset < 0
        start    = find((abs(StreamData.rpm{i}(2:end-1)) < est_avg_rpm - Trig),1);          %first index at which trig occurs
    else
        start    = find((abs(StreamData.rpm{i}(2:end-1)) > est_avg_rpm + Trig),1);          %first index at which trig occurs
    end
        
%     while start == 1
%         StreamData.rpm{i} = fcleanup(StreamData.rpm{i}, 'smoothdata', 'loess', 700); 
%         if offset < 0
%             start    = find((abs(StreamData.rpm{i}(2:end)) < est_avg_torque - Trig),1);          %first index at which trig occurs
%         else
%             start    = find((abs(StreamData.rpm{i}(2:end)) > est_avg_torque + Trig),1);          %first index at which trig occurs
%         end
%     end
    
    start = start - 500; %adjust for delay in rpm trig due to smoothing
    
    rpm_motor0  = mean(StreamData.rpm{i}(1:start));                         %avg motor rpm prior to trig
    dazds       = rpm_motor0/60*360/SR;                                     %rate of change in az over samples
    az_trig     = dazds*start;                                              %estimated az at the trig location [deg]
    
    %create ref_angle with slew
    samples = SR*offset/slew; %samples to complete slew
    comm_angle{i}   = dazds*(0:length(StreamData.encoder{i})-1);
    ref_angle{i}   = dazds*(0:length(StreamData.encoder{i})-1);
    for idx_ref = start:(start+samples)
        comm_angle{i}(idx_ref) = comm_angle{i}(idx_ref) + slew/SR*(idx_ref-start);
    end
    comm_angle{i}(start+samples+1:end) = comm_angle{i}(start+samples+1:end) + offset;
    ref_angle{i}   =  StreamData.unwrap{i} - ref_angle{i};


    %calc angle err
    idx = 1;
    while idx <= length(comm_angle{i})
        StreamData.angle_err{i}(idx) = comm_angle{i}(idx) - StreamData.unwrap{i}(idx);
        idx = idx+1;
    end
    reset_az = StreamData.unwrap{i} - StreamData.unwrap{i}(start); 
    
    if (start - pre_trig*SR) < 0; idx_trig = 1:(start + post_trig*SR);
    else; idx_trig = (start - pre_trig*SR):(start + post_trig*SR); %index .25s pre-trig and 1s post-trig
    end

    time            = Time - Time(start); time = time(idx_trig);%get time vect
    time0           = Time - Time(start); time0 = time0(idx_trig - 1);
    Mz_trans(:,i)   = StreamData.Mz_inner{i}(idx_trig);
    Q_est_trans(:,i)   = fcleanup(StreamData.IQ{i}(idx_trig)*Kt/sqrt(2),'smoothdata','loess',500); %estimated motor torque based on IQ and Kt
    RPM_dazdt(:,i)  = StreamData.rpm{i}(idx_trig);
    Angle_err(:,i)  = StreamData.angle_err{i}(idx_trig);
    Ref_angle(:,i)     = ref_angle{i}(idx_trig);
    Meas_angle(:,i) = reset_az(idx_trig);
    
    T           = fcleanup(StreamData.Fz_inner{i}, 'smoothdata', 'loess', 700);
    T_trans(:,i)  = T(idx_trig);
    
    PhaseSync.Meas_angle(:,i) = Meas_angle(:,i);
    PhaseSync.time = time; %get time vect
    PhaseSync.Thrust(:,i) = T_trans(:,i);
    PhaseSync.Angle_err(:,i) = StreamData.angle_err{i}(idx_trig);
    PhaseSync.Ref_angle(:,i) = Ref_angle(:,i);
    PhaseSync.Speed(:,i) = RPM_dazdt(:,i);
    PhaseSync.Torque(:,i) = Mz_trans(:,i)/GR;
    PhaseSync.Q_est(:,i) = Q_est_trans(:,i);
end

if size(PhaseSync.Thrust,2) == 1
    PhaseSync.rev = PhaseSync.Meas_angle';
    PhaseSync.T_avg = PhaseSync.Thrust';
    PhaseSync.T_err = tinv(.975,1) * Fb * ones(1,length(PhaseSync.Thrust));
    PhaseSync.Ang_err_avg = PhaseSync.Angle_err';
    PhaseSync.Ang_err_err = tinv(.975,1) * Az_b * ones(1,length(PhaseSync.Angle_err));
    PhaseSync.ref_ang_avg = PhaseSync.Ref_angle';
    PhaseSync.ref_ang_err = tinv(.975,1) * Az_b * ones(1,length(PhaseSync.Ref_angle));
    PhaseSync.Speed_avg = PhaseSync.Speed';
    PhaseSync.Speed_err = 0 * ones(1,length(PhaseSync.Speed));
    PhaseSync.Torque_avg = PhaseSync.Torque';
    PhaseSync.Torque_err = tinv(.975,1) * Mb * ones(1,length(PhaseSync.Torque));
    PhaseSync.Q_est_avg = PhaseSync.Q_est';
    PhaseSync.Q_est_err = tinv(.975,1) * (IQ_b*Kt/sqrt(2)) * ones(1,length(PhaseSync.Q_est));
        
elseif size(PhaseSync.Thrust,2) > 1
    PhaseSync.rev = mean(PhaseSync.Meas_angle');
    PhaseSync.T_avg = mean(PhaseSync.Thrust');
    PhaseSync.T_err = tinv(.975,size(PhaseSync.Thrust,2)) * sqrt(std(PhaseSync.Thrust').^2 + Fb^2)/sqrt(size(PhaseSync.Thrust,2));
    PhaseSync.Ang_err_avg = mean(PhaseSync.Angle_err');
    PhaseSync.Ang_err_err = tinv(.975,size(PhaseSync.Angle_err,2)) * sqrt(std(PhaseSync.Angle_err').^2 + Az_b^2)/sqrt(size(PhaseSync.Angle_err,2));
    PhaseSync.ref_ang_avg = mean(PhaseSync.Ref_angle');
    PhaseSync.ref_ang_err = tinv(.975,size(PhaseSync.Ref_angle,2)) * sqrt(std(PhaseSync.Ref_angle').^2 + Az_b^2)/sqrt(size(PhaseSync.Ref_angle,2));
    PhaseSync.Speed_avg = mean(PhaseSync.Speed');
    PhaseSync.Speed_err = tinv(.975,size(PhaseSync.Speed,2)) * std(PhaseSync.Speed')/sqrt(size(PhaseSync.Speed,2));
    PhaseSync.Torque_avg = mean(PhaseSync.Torque');
    PhaseSync.Torque_err = tinv(.975,size(PhaseSync.Torque,2)) * sqrt(std(PhaseSync.Torque').^2 + Mb^2)/sqrt(size(PhaseSync.Torque,2));
    PhaseSync.Q_est_avg = mean(PhaseSync.Q_est');
    PhaseSync.Q_est_err = tinv(.975,size(PhaseSync.Q_est,2)) * sqrt(std(PhaseSync.Q_est').^2 + (IQ_b*Kt/sqrt(2))^2)/sqrt(size(PhaseSync.Q_est,2));
end

fprintf('\nPhase-sync processing done.\n')

end

