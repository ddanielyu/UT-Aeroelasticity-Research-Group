function [PhaseSync] = runPhaseSync_dual(StreamData,phaseSync_test,offset,Trig,GR,slew)
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

Fs = 4e3; %Sample rate [Hz]
pre_trig = 0.25;
post_trig = 1;

%% Find New Stream Data
i = 1; %counter for indexing each test within a set

for i = 1:length(phaseSync_test)
    fprintf('\n%s%s\n','Processing',phaseSync_test{i})        

    
    Trig = 4*std(StreamData.Mz_inner{i});
    
    StreamData.Mz_inner{i} = StreamData.Mz_inner{i};
    
    Time     = 0:1/Fs:(length(StreamData.encoder{i})-1)/Fs;
    
    %find angle trigger and calculate the rev-averaged angle error
    est_avg_torque = abs(nanmean(StreamData.Mz_inner{i}));
    
    if offset < 0
        start    = find((abs(StreamData.Mz_inner{i}(2:end-1)) < est_avg_torque - Trig),1);          %first index at which trig occurs
    else
        start    = find((abs(StreamData.Mz_inner{i}(2:end-1)) > est_avg_torque + Trig),1);          %first index at which trig occurs
    end
        
%     while start == 1
%         StreamData.rpm{idx_PS} = fcleanup(StreamData.rpm{idx_PS}, 'smoothdata', 'loess', 700); 
%         if offset < 0
%             start    = find((abs(StreamData.rpm{idx_PS}(2:end)) < est_avg_torque - Trig),1);          %first index at which trig occurs
%         else
%             start    = find((abs(StreamData.rpm{idx_PS}(2:end)) > est_avg_torque + Trig),1);          %first index at which trig occurs
%         end
%     end
    
%     start = start; %adjust for delay in rpm trig due to smoothing

    %Filter Data 
%     StreamData.rpm{i} = fcleanup(StreamData.rpm{i},'smoothdata', 'loess', 700);    
%     StreamData.accel{i} = fcleanup(StreamData.accel{i},'smoothdata', 'loess', 700);
%     StreamData.Mz_inner{i} = fcleanup(StreamData.Mz_inner{i},'smoothdata','loess',250);
%     StreamData.Fz_inner{i} = fcleanup(StreamData.Fz_inner{i}, 'smoothdata', 'loess', 700);
%     StreamData.IQ{i} = fcleanup(StreamData.IQ{i},'smoothdata','loess',500); %estimated motor torque based on IQ and Kt
    StreamData.rpm{i} = savitzkyGolayFilt(StreamData.rpm{i},2,0,401);
    StreamData.accel{i} = savitzkyGolayFilt(StreamData.accel{i},2,0,401);
    StreamData.accel{i} = savitzkyGolayFilt(StreamData.accel{i},2,0,401);
    StreamData.Mz_inner{i} = savitzkyGolayFilt(StreamData.Mz_inner{i},2,0,401);
    dTdt        = (StreamData.Fz_inner{i}(2:end) - StreamData.Fz_inner{i}(1:end-1))*Fs;
    StreamData.Fz_inner{i} = savitzkyGolayFilt(StreamData.Fz_inner{i},2,0,401);
    StreamData.IQ{i} = savitzkyGolayFilt(StreamData.IQ{i},2,0,401); %estimated motor torque based on IQ and Kt
    
    rpm_motor0  = mean(StreamData.rpm{i}(1:start));                         %avg motor rpm prior to trig

    ct_den = StreamData.rho{i} * (pi * StreamData.R^2) * (rpm_motor0*2*pi/60*StreamData.R).^2;
    cts_bias = 18.75 / ct_den / StreamData.sigma;

    dazds       = rpm_motor0/60*360/Fs;                                     %rate of change in az over samples
    az_trig     = dazds*start;                                              %estimated az at the trig location [deg]
    
    %create ref_angle with slew
    samples = Fs*offset/slew; %samples to complete slew
    comm_angle{i}   = dazds*(0:length(StreamData.encoder{i})-1);
    ref_angle   = dazds*(0:length(StreamData.encoder{i})-1);
    for idx_ref = start:(start+samples)
        comm_angle{i}(idx_ref) = comm_angle{i}(idx_ref) + slew/Fs*(idx_ref-start);
    end
    comm_angle{i}(start+samples+1:end) = comm_angle{i}(start+samples+1:end) + offset;
    ref_angle   =  StreamData.unwrap{i} - ref_angle' - StreamData.unwrap{i}(1);


    %calc angle err
    idx = 1;
    while idx <= length(comm_angle{i})
        StreamData.angle_err{i}(idx) = comm_angle{i}(idx) - StreamData.unwrap{i}(idx);
        idx = idx+1;
    end
    hi = StreamData.unwrap{i} - StreamData.unwrap{i}(start); 
    
    if (start - pre_trig*Fs) < 0; idx_trig = 1:(start + post_trig*Fs);
    else; idx_trig = (start - pre_trig*Fs):(start + post_trig*Fs); %index .25s pre-trig and 1s post-trig
    end
    
    if start < Fs/4; continue; end
    

    time            = Time - Time(start); time = time(idx_trig);%get time vect
    time0           = Time - Time(start); time0 = time0(idx_trig - 1);
    Mz_trans(:,i)   = StreamData.Mz_inner{i}(idx_trig);
    Q_est_trans(:,i)   = StreamData.IQ{i}(idx_trig)*Kt/sqrt(2); %estimated motor torque based on IQ and Kt
    RPM_dazdt(:,i)  = StreamData.rpm{i}(idx_trig);
    Angle_err(:,i)  = StreamData.angle_err{i}(idx_trig);
    Ref_angle(:,i)     = ref_angle(idx_trig);
    Meas_angle(:,i) = hi(idx_trig);
    
    T           = StreamData.Fz_inner{i};
    T_trans     = T(idx_trig);
    
    %calc dT/dt to compare with Thrust variations due to RPM
    dTdt        = savitzkyGolayFilt(dTdt,2,0,401);
    dTdt        = savitzkyGolayFilt(dTdt,2,0,401);
    dTdt_trans  = dTdt(idx_trig);
    
    PhaseSync.Meas_angle(:,i) = Meas_angle(:,i); %get time vec
    PhaseSync.time = time; %get time vect
    PhaseSync.Thrust(:,i) = T_trans;
    PhaseSync.dTdt(:,i) = dTdt_trans;
    PhaseSync.Angle_err(:,i) = StreamData.angle_err{i}(idx_trig);
    PhaseSync.Ref_angle(:,i) = Ref_angle(:,i);
    PhaseSync.Speed(:,i) = RPM_dazdt(:,i);
    PhaseSync.Torque(:,i) = Mz_trans(:,i)/GR;
    PhaseSync.Q_est(:,i) = Q_est_trans(:,i);
    PhaseSync.Cts(:,i) = T_trans./(StreamData.rho{i}*pi*StreamData.R^2.*(PhaseSync.Speed(:,i)*2*pi/60*StreamData.R).^2*StreamData.sigma);
    PhaseSync.accel(:,i) = StreamData.accel{i}(idx_trig);
end

PhaseSync.rev = mean(PhaseSync.Meas_angle');
PhaseSync.T_avg = mean(PhaseSync.Thrust');
PhaseSync.T_err = tinv(.975,size(PhaseSync.Thrust,2)) * sqrt(std(PhaseSync.Thrust').^2 + Fb^2)/sqrt(size(PhaseSync.Thrust,2));
PhaseSync.dTdt_avg = mean(PhaseSync.dTdt');
PhaseSync.dTdt_err = tinv(.975,size(PhaseSync.dTdt,2)) * std(PhaseSync.dTdt')/sqrt(size(PhaseSync.dTdt,2));
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
PhaseSync.Cts_avg = mean(PhaseSync.Cts');
PhaseSync.Cts_err = tinv(.975,size(PhaseSync.Cts,2)) * sqrt(std(PhaseSync.Cts').^2 + cts_bias^2)/sqrt(size(PhaseSync.Cts,2));
PhaseSync.accel_avg = mean(PhaseSync.accel');
PhaseSync.accel_err = tinv(.975,size(PhaseSync.accel,2)) * std(PhaseSync.accel')/sqrt(size(PhaseSync.accel,2));

fprintf('\nPhase-sync processing done.\n')

end

