function [PhaseSync] = runPhaseSync(StreamData,phaseSync_test,offset,GR,slew,motor)
%{

This function processes data from phase-sync tests from Dual Motor testing
2021-2022 Testing and outputs all phase-sync data into PhaseSync struc

WRITTEN BY: MATT ASPER
DATE: MARCH 17, 2022

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

Fs = 4e3; %Sample rate [Hz]
pre_trig = 0.25;
post_trig = 1;

%% Find New Stream Data
ij = 1; %counter for indexing each test within a set

for i = 1:length(phaseSync_test)
    current_test_name = split(phaseSync_test{i},{'.'});
    fprintf('\n%s%s\n','Processing ',current_test_name{1}) 
    
    %Find phase-sync file in StreamData

    Trig = 1; %1 deg change in index
        
    Time     = 0:1/Fs:(length(StreamData.encoder1{i})-1)/Fs;
    
    %find angle trigger and calculate the rev-averaged angle error
    avg_index = nanmean(StreamData.index{i}(1:Fs/2)); %avg first 0.5s of index data
    
    if offset > 0
        start    = find(StreamData.index{i} > avg_index + Trig,1);          %first index at which trig occurs
    elseif offset < 0
        start    = find(StreamData.index{i} < avg_index - Trig,1);          %first index at which trig occurs
    end
    
    start = start - Trig/slew*Fs; %re-align data with correct 0
    start = int64(start);
    
    RPM1 = (StreamData.encoder1{i}(2:end) - StreamData.encoder1{i}(1:end-1))*Fs/360*60;
    RPM1 = fcleanup(RPM1, 'smoothdata', 'loess', 700);
    
    RPM2 = (StreamData.encoder2{i}(2:end) - StreamData.encoder2{i}(1:end-1))*Fs/360*60;
    RPM2 = fcleanup(RPM2, 'smoothdata', 'loess', 700);
    
    %avg motor rpm prior to trig
    rpm_motor = mean(RPM1(1:start));
    rpm_motor2 = mean(RPM2(1:start));

    dazds       = rpm_motor/60*360/Fs;                                     %rate of change in az over samples
    dazds2       = rpm_motor2/60*360/Fs;
    
    %create ref_angle with slew
    samples = Fs*offset/slew; %samples to complete slew
    
    ref_angle   = dazds*(0:length(StreamData.encoder1{i})-1); %ref angle for servo motor
    ref_angle2   = dazds2*(0:length(StreamData.encoder2{i})-1); %ref angle for follower motor
    
    %commanded angle (with slew)
    if motor == 'f'
        comm_angle{i}   = dazds2*(0:length(StreamData.encoder1{i})-1);
    elseif motor == 's'
        comm_angle{i}   = dazds*(0:length(StreamData.encoder1{i})-1);
    end
    for idx_ref = start:(start+samples)
        comm_angle{i}(idx_ref) = comm_angle{i}(idx_ref) + slew/Fs*(idx_ref-start);
    end
    comm_angle{i}(start+samples+1:end) = comm_angle{i}(start+samples+1:end) + offset;
    
    
    ref_angle   =  StreamData.encoder1{i} - ref_angle'; %servo
    ref_angle2   =  StreamData.encoder2{i} - ref_angle2'; %follower


    %calc angle err
    idx = 1;
    while idx <= length(comm_angle{i})
        if motor == 'f'
            StreamData.angle_err1{i}(idx) = ref_angle(idx) - StreamData.encoder1{i}(idx);
            StreamData.angle_err2{i}(idx) = comm_angle{i}(idx) - StreamData.encoder2{i}(idx);
        elseif motor == 's'
            StreamData.angle_err1{i}(idx) = comm_angle{i}(idx) - StreamData.encoder1{i}(idx);
            StreamData.angle_err2{i}(idx) = ref_angle2(idx) - StreamData.encoder2{i}(idx);
        end    
        idx = idx+1;
    end
    bias_angle1 = StreamData.encoder1{i} - StreamData.encoder1{i}(start); 
    bias_angle2 = StreamData.encoder2{i} - StreamData.encoder2{i}(start); 
    
    idx_trig = (start - pre_trig*Fs):(start + post_trig*Fs); %index .25s pre-trig and 1s post-trig
    

    time            = Time - Time(start); time = time(idx_trig);%get time vect
    time0           = Time - Time(start); time0 = time0(idx_trig - 1);
    Mz_outer_trans(:,i)   = StreamData.Mz_outer{i}(idx_trig);
    Mz_inner_trans(:,i)   = StreamData.Mz_inner{i}(idx_trig);
    Q1_est_trans(:,i)   = fcleanup(StreamData.IQ1{i}(idx_trig)*Kt/sqrt(2)*GR,'smoothdata','loess',300); %estimated servo rotor torque
    Q2_est_trans(:,i)   = fcleanup(StreamData.IQ2{i}(idx_trig)*Kt/sqrt(2)*GR,'smoothdata','loess',300); %estimated follower rotor torque
    RPM1_trans            = RPM1(idx_trig);
    RPM2_trans            = RPM2(idx_trig);
    Angle_err1(:,i)  = StreamData.angle_err1{i}(idx_trig);
    Angle_err2(:,i)  = StreamData.angle_err2{i}(idx_trig);
    Ref_angle1(:,i)     = ref_angle(idx_trig);
    Ref_angle2(:,i)     = ref_angle2(idx_trig);
    Meas_angle1(:,i) = bias_angle1(idx_trig);
    Meas_angle2(:,i) = bias_angle2(idx_trig);

    T_outer_trans(:,i)   = abs(StreamData.Fz_outer{i}(idx_trig));
    T_outer_trans(:,i)  = fcleanup(T_outer_trans(:,i), 'smoothdata', 'loess', 500);

    T_inner_trans(:,i)   = abs(StreamData.Fz_inner{i}(idx_trig));
    T_inner_trans(:,i)  = fcleanup(T_inner_trans(:,i), 'smoothdata', 'loess', 500);

    PhaseSync.Meas_angle1(:,ij) = Meas_angle1(:,i); %get time vec
    PhaseSync.Meas_angle2(:,ij) = Meas_angle2(:,i); %get time vec
    PhaseSync.time = time; %get time vect
    PhaseSync.Thrust_outer(:,ij) = T_outer_trans(:,i);
    PhaseSync.Thrust_inner(:,ij) = T_inner_trans(:,i);
    PhaseSync.Angle_err1(:,ij) = StreamData.angle_err1{i}(idx_trig);
    PhaseSync.Angle_err2(:,ij) = StreamData.angle_err2{i}(idx_trig);
    PhaseSync.Ref_angle1(:,ij) = Ref_angle1(:,i);
    PhaseSync.Ref_angle2(:,ij) = Ref_angle2(:,i);
    PhaseSync.index(:,ij) = StreamData.index{i}(idx_trig);
    PhaseSync.servo_speed(:,ij) = RPM1_trans;
    PhaseSync.follower_speed(:,ij) = RPM2_trans;
    PhaseSync.Torque_outer(:,ij) = Mz_outer_trans(:,i);
    PhaseSync.Torque_inner(:,ij) = Mz_inner_trans(:,i);
    PhaseSync.Q1_est(:,ij) = Q1_est_trans(:,i);
    PhaseSync.Q2_est(:,ij) = Q2_est_trans(:,i);
    ij = ij + 1;
end

PhaseSync.rev1 = mean(PhaseSync.Meas_angle1');
PhaseSync.rev2 = mean(PhaseSync.Meas_angle2');
PhaseSync.T_outer_avg = mean(PhaseSync.Thrust_outer');
PhaseSync.T_outer_err = tinv(.975,size(PhaseSync.Thrust_outer,2)) * sqrt(std(PhaseSync.Thrust_outer').^2 + Fb^2)/sqrt(size(PhaseSync.Thrust_outer,2));
PhaseSync.T_inner_avg = mean(PhaseSync.Thrust_inner');
PhaseSync.T_inner_err = tinv(.975,size(PhaseSync.Thrust_inner,2)) * sqrt(std(PhaseSync.Thrust_inner').^2 + Fb^2)/sqrt(size(PhaseSync.Thrust_inner,2));
PhaseSync.Ang_err1_avg = mean(PhaseSync.Angle_err1');
PhaseSync.Ang_err1_err = tinv(.975,size(PhaseSync.Angle_err1,2)) * sqrt(std(PhaseSync.Angle_err1').^2 + Az_b^2)/sqrt(size(PhaseSync.Angle_err1,2));
PhaseSync.Ang_err2_avg = mean(PhaseSync.Angle_err2');
PhaseSync.Ang_err2_err = tinv(.975,size(PhaseSync.Angle_err2,2)) * sqrt(std(PhaseSync.Angle_err2').^2 + Az_b^2)/sqrt(size(PhaseSync.Angle_err2,2));
PhaseSync.ref_ang1_avg = mean(PhaseSync.Ref_angle1');
PhaseSync.ref_ang1_err = tinv(.975,size(PhaseSync.Ref_angle1,2)) * sqrt(std(PhaseSync.Ref_angle1').^2 + Az_b^2)/sqrt(size(PhaseSync.Ref_angle1,2));
PhaseSync.ref_ang2_avg = mean(PhaseSync.Ref_angle2');
PhaseSync.ref_ang2_err = tinv(.975,size(PhaseSync.Ref_angle2,2)) * sqrt(std(PhaseSync.Ref_angle2').^2 + Az_b^2)/sqrt(size(PhaseSync.Ref_angle2,2));
PhaseSync.servo_speed_avg = mean(PhaseSync.servo_speed');
PhaseSync.servo_speed_err = tinv(.975,size(PhaseSync.servo_speed,2)) * std(PhaseSync.servo_speed')/sqrt(size(PhaseSync.servo_speed,2));
PhaseSync.follower_speed_avg = mean(PhaseSync.follower_speed');
PhaseSync.follower_speed_err = tinv(.975,size(PhaseSync.follower_speed,2)) * std(PhaseSync.follower_speed')/sqrt(size(PhaseSync.follower_speed,2));
PhaseSync.Torque_outer_avg = mean(PhaseSync.Torque_outer');
PhaseSync.Torque_outer_err = tinv(.975,size(PhaseSync.Torque_outer,2)) * sqrt(std(PhaseSync.Torque_outer').^2 + Mb^2)/sqrt(size(PhaseSync.Torque_outer,2));
PhaseSync.Torque_inner_avg = mean(PhaseSync.Torque_inner');
PhaseSync.Torque_inner_err = tinv(.975,size(PhaseSync.Torque_inner,2)) * sqrt(std(PhaseSync.Torque_inner').^2 + Mb^2)/sqrt(size(PhaseSync.Torque_inner,2));
PhaseSync.Q1_est_avg = mean(PhaseSync.Q1_est');
PhaseSync.Q1_est_err = tinv(.975,size(PhaseSync.Q1_est,2)) * sqrt(std(PhaseSync.Q1_est').^2 + (IQ_b*Kt/sqrt(2)*GR)^2)/sqrt(size(PhaseSync.Q1_est,2));
PhaseSync.Q2_est_avg = mean(PhaseSync.Q2_est');
PhaseSync.Q2_est_err = tinv(.975,size(PhaseSync.Q2_est,2)) * sqrt(std(PhaseSync.Q2_est').^2 + (IQ_b*Kt/sqrt(2)*GR)^2)/sqrt(size(PhaseSync.Q2_est,2));
PhaseSync.index_avg = mean(PhaseSync.index');
PhaseSync.index_err = tinv(.975,size(PhaseSync.index,2)) * sqrt(std(PhaseSync.index').^2 + Az_b^2)/sqrt(size(PhaseSync.index,2));

fprintf('\nPhase-sync processing done.\n')

end

