%{
This script processes data from Aug 8, 2022-Present Tripod Stand Testing. Either
steady-state, phase-sync, or dynamic RPM data can be processed.

This processing code uses the 220212 Labview code (single motor)

Written By: Matt Asper
Date: 20 Sept 2022

%}

clc; clear; close all;
load('colors.mat');

%% Inputs
files_dir = uigetdir(); %directory of .csv data files
loads = input('Plot with loads? [y/n]: ','s');
source_dir = pwd; %directory of MATLAB scripts
SR = input('Sampling Frequency [Hz]: ');

%% Constants
%smoothing parameters
rpm_smoothing = 700; %samples to smooth rpm signal
torque_smoothing = 1500; %samples to smooth ESTIMATED torque signal
thrust_smoothing = 700; %samples to smooth thrust signal
accel_smoothing = 700; %samples to smooth acceleration signal after differentiating rpm

%% Load Data
conditions = [54	29.88]; % [T(Farenh), % humidity, P(in.Hg)]

[mdata,MeanData,steady_test,phaseSync_test,dynRPM_test] = loadFiles(source_dir,files_dir,conditions); %load files

%Extract variables
collective = mdata.MeanCollective(1);
% collective = input('Collective [deg]: ');
GR = mdata.GearRatio(1);
% GR = input('Gear Ratio: ');


[MeanData,StreamData] = loadStreamTripod(mdata,MeanData,source_dir,files_dir,collective,SR,GR,rpm_smoothing); %load streamdata

%% Process ALL Data
fprintf('\nProcessing data sets...\n')

%Execute remainder of LoadsProcessing.m code
[StreamData,SortedData] = fSortStream(StreamData, SR, GR);

RevData = fRevolutionAvg(SortedData);
AvgData = fTotalAvg(RevData,SortedData,StreamData);


%% Process Steady Data
%clear loads to avoid corrupting plots
if loads == 'n'
    for i = 1:length(AvgData.avg_cts_inner)
        AvgData.avg_cts_inner{i} = NaN;
        AvgData.err_cts_inner{i} = NaN;
        AvgData.avg_cps_inner{i} = NaN;
        AvgData.err_cps_inner{i} = NaN;
    end
end

fprintf('\nProcessing steady-only data...\n')

if isempty(steady_test) == 0
    [Averages] = runSteady(StreamData,AvgData,steady_test);
end

%% Process Phase Sync

Trig = 10; %rpm spike from nominal rpm to align data sets for phase-sync

fprintf('\nRunning phase-sync processing...\n')
if isempty(phaseSync_test) == 0
    offset = input('Offset Angle [deg]: ');
    slew = input('Angle Slew Rate [deg/s]: ');
    [PhaseSync] = runPhaseSync(StreamData,phaseSync_test,offset,Trig,GR,slew,SR);

    %clear loads to avoid corrupting plots
    if loads == 'n'
        for i = 1:length(PhaseSync.T_avg)
            PhaseSync.T_avg = NaN;
            PhaseSync.T_err = NaN;
            PhaseSync.Torque_avg = NaN;
            PhaseSync.Torque_err = NaN;
        end
    end
    
end


%% Process Dynamic RPM

Trig = 50; %rpm change from nominal rpm to align data sets for dynamic-rpm testing

fprintf('\nRunning dynamic RPM processing...\n')
if isempty(dynRPM_test) == 0
    slew = input('Ramp Rate [RPM/s]: ');
    [DynRPM] = runDynRPM(StreamData,dynRPM_test,Trig,GR,slew,SR,torque_smoothing,thrust_smoothing,accel_smoothing,rpm_smoothing);

    %clear loads to avoid corrupting plots
    if loads == 'n'
        for i = 1:length(DynRPM.T_avg)
            DynRPM.T_avg = NaN;
            DynRPM.T_err = NaN;
            DynRPM.Torque_avg = NaN;
            DynRPM.Torque_err = NaN;
        end
    end
    
end
%% Plotting
close all; clc;

if isempty(steady_test) == 0; [f1,f2,f3,f4,f5] = plotSteady(Averages,collective); end
if isempty(phaseSync_test) == 0; [f6,f7] = plotPhaseSync(PhaseSync,loads); end
if isempty(dynRPM_test) == 0; [f8,f9,f10,f11] = plotDynRPM(DynRPM,loads,GR); end


%% Saving
save_dir = uigetdir();
if isempty(steady_test) == 0
    save(fullfile(save_dir,'Data'),'Averages');
    saveas(f1,fullfile(save_dir,f1.Name),'jpg')
    saveas(f2,fullfile(save_dir,f2.Name),'jpg')
    saveas(f3,fullfile(save_dir,f3.Name),'jpg')
    saveas(f4,fullfile(save_dir,f4.Name),'jpg')
    saveas(f5,fullfile(save_dir,f5.Name),'jpg')
elseif isempty(phaseSync_test) == 0
    save(fullfile(save_dir,'Data'),'PhaseSync');
    saveas(f6,fullfile(save_dir,f6.Name),'jpg')
    saveas(f7,fullfile(save_dir,f7.Name),'jpg')
elseif isempty(dynRPM_test) == 0
    save(fullfile(save_dir,'Data'),'DynRPM');
    saveas(f8,fullfile(save_dir,f8.Name),'jpg')
    saveas(f9,fullfile(save_dir,f9.Name),'jpg')
    saveas(f10,fullfile(save_dir,f10.Name),'jpg')
    saveas(f11,fullfile(save_dir,f11.Name),'jpg')
end


