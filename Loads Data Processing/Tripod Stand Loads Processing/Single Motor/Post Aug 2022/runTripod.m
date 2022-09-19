%{
This script processes data from Aug 8, 2022-Present Tripod Stand Testing. Either
steady-state or phase-sync data can be processed.

This processing code uses the 220212 Labview code (single motor)

Written By: Matt Asper
Date: 08 Aug 2022

%}

clc; clear; close all;
load('colors.mat');

%% Constants

Trig = 10; %rpm spike from nominal rpm to find phase sync trigger

%% Inputs
files_dir = uigetdir(); %directory of .csv data files
loads = input('Plot with loads? [y/n]: ','s');
source_dir = pwd; %directory of MATLAB scripts
SR = input('Sampling Frequency [Hz]: ');

%% Load Data
conditions = [54	29.88]; % [T(Farenh), % humidity, P(in.Hg)]

[mdata,MeanData,steady_test,phaseSync_test] = loadFiles(source_dir,files_dir,conditions); %load files

%Extract variables
% collective = mdata.MeanCollective(1);
collective = input('Collective [deg]: ');
% GR = mdata.GearRatio(1);
GR = input('Gear Ratio: ');


[MeanData,StreamData] = loadStreamTripod(mdata,MeanData,source_dir,files_dir,collective,SR,GR); %load streamdata

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

%% Plotting
close all; clc;

if isempty(steady_test) == 0; [f1,f2,f3,f4] = plotSteady(Averages,collective); end
if isempty(phaseSync_test) == 0; [f5,f6] = plotPhaseSync(PhaseSync,loads); end


%% Saving
save_dir = uigetdir();
if isempty(steady_test) == 0
    save(fullfile(save_dir,'Data'),'Averages');
    saveas(f1,fullfile(save_dir,f1.Name),'jpg')
    saveas(f2,fullfile(save_dir,f2.Name),'jpg')
    saveas(f3,fullfile(save_dir,f3.Name),'jpg')
    saveas(f4,fullfile(save_dir,f4.Name),'jpg')
elseif isempty(phaseSync_test) == 0
    save(fullfile(save_dir,'Data'),'PhaseSync');
    saveas(f5,fullfile(save_dir,f5.Name),'jpg')
    saveas(f6,fullfile(save_dir,f6.Name),'jpg')
end


