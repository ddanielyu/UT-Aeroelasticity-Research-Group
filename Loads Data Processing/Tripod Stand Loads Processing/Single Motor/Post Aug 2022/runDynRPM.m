function [DynRPM] = runDynRPM(StreamData,dynRPM_test,Trig,GR,slew,SR,torque_smoothing,thrust_smoothing)
%{

This function processes data from dynamic RPM tests from Tripod Stand
Post Aug-2022 Testing and outputs all data into DynRPM struc

DynRPM:
time
rev (motor revolution!!)
Thrust (sets separated by columns); T_avg; T_err;
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
post_trig = 2;

%% Find New Stream Data

for i = 1:length(dynRPM_test)

    Time     = 0:1/SR:(length(StreamData.encoder{i})-1)/SR;
    
    %estimate the avg rpm based on first .5s of data
    est_avg_rpm = abs(nanmean(StreamData.rpm{i}(1:SR/2)));
    
    %first index at which trig occurs
    start    = find((abs(StreamData.rpm{i}(2:end-1)) > est_avg_rpm + Trig),1);
    
    dt  = Trig/slew; %time for motor to manuever to "Trig" RPM based on slew rate
    start = start - dt*SR - 700; %adjust start location due to slew rate AND smoothing

    %create indexed-array and adjust indicies if Pre- and Post- Trig durations are not appropriate
    if (start - pre_trig*SR) < 0; idx_trig = 1:(start + post_trig*SR);
    else; idx_trig = (start - pre_trig*SR):(start + post_trig*SR);
    end
    
    reset_az = StreamData.unwrap{i} - StreamData.unwrap{i}(start); 
    
    
    %-----------Transient Data--------------%
    time                = Time - Time(start); time = time(idx_trig);%get time vect
    time0               = Time - Time(start); time0 = time0(idx_trig - 1);
    
    %Measured torque
    Mz_trans(:,i)       = StreamData.Mz_inner{i}(idx_trig);
    
    %estimated motor torque based on IQ and Kt
    Q_est_trans(:,i)    = fcleanup(StreamData.IQ{i}(idx_trig)*Kt/sqrt(2),'smoothdata','loess',torque_smoothing);
    
    %transient speed
    RPM_dazdt(:,i)      = StreamData.rpm{i}(idx_trig);
    
    %motor revolution
    Meas_angle(:,i)     = reset_az(idx_trig);
    
    %thrust
    T                   = fcleanup(StreamData.Fz_inner{i}, 'smoothdata', 'loess', thrust_smoothing);
    T_trans(:,i)        = T(idx_trig);
    
    %compile data
    DynRPM.Meas_angle(:,i) = Meas_angle(:,i);
    DynRPM.time = time; %get time vect
    DynRPM.Thrust(:,i) = T_trans(:,i);
    DynRPM.Speed(:,i) = RPM_dazdt(:,i);
    DynRPM.Torque(:,i) = Mz_trans(:,i)/GR;
    DynRPM.Q_est(:,i) = Q_est_trans(:,i);
end

%Average data sets and calculate error with 95% CI
if size(DynRPM.Thrust,2) == 1
    DynRPM.rev = DynRPM.Meas_angle';
    DynRPM.T_avg = DynRPM.Thrust';
    DynRPM.T_err = tinv(.975,1) * Fb * ones(1,length(DynRPM.Thrust));
    DynRPM.Speed_avg = DynRPM.Speed';
    DynRPM.Speed_err = 0 * ones(1,length(DynRPM.Speed));
    DynRPM.Torque_avg = DynRPM.Torque';
    DynRPM.Torque_err = tinv(.975,1) * Mb * ones(1,length(DynRPM.Torque));
    DynRPM.Q_est_avg = DynRPM.Q_est';
    DynRPM.Q_est_err = tinv(.975,1) * (IQ_b*Kt/sqrt(2)) * ones(1,length(DynRPM.Q_est));
        
elseif size(DynRPM.Thrust,2) > 1
    DynRPM.rev = mean(DynRPM.Meas_angle');
    DynRPM.T_avg = mean(DynRPM.Thrust');
    DynRPM.T_err = tinv(.975,size(DynRPM.Thrust,2)) * sqrt(std(DynRPM.Thrust').^2 + Fb^2)/sqrt(size(DynRPM.Thrust,2));
    DynRPM.Speed_avg = mean(DynRPM.Speed');
    DynRPM.Speed_err = tinv(.975,size(DynRPM.Speed,2)) * std(DynRPM.Speed')/sqrt(size(DynRPM.Speed,2));
    DynRPM.Torque_avg = mean(DynRPM.Torque');
    DynRPM.Torque_err = tinv(.975,size(DynRPM.Torque,2)) * sqrt(std(DynRPM.Torque').^2 + Mb^2)/sqrt(size(DynRPM.Torque,2));
    DynRPM.Q_est_avg = mean(DynRPM.Q_est');
    DynRPM.Q_est_err = tinv(.975,size(DynRPM.Q_est,2)) * sqrt(std(DynRPM.Q_est').^2 + (IQ_b*Kt/sqrt(2))^2)/sqrt(size(DynRPM.Q_est,2));
end

fprintf('\nDynamic RPM processing done.\n')

end

