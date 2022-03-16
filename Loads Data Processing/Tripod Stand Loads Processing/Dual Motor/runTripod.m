%{
This script processes data from Fall 2021 Tripod Stand Testing. Either
steady-state or phase-sync data can be processed.

This processing code is a modification to the "LoadsProcessing.m" code, and
essentially uses only fLoadData

Written By: Matt Asper
Date: 17 Feb 2022

%}

clc; clear; close all;
load('colors.mat');

%% Constants

Trig = 6; %rpm spike from nominal rpm to find phase sync trigger
source_dir = pwd; %directory of MATLAB scripts

%% Inputs
files_dir = uigetdir(); %directory of .csv data files
loads = input('Plot with loads? [y/n]: ','s');

%% Load Data
conditions = [54	29.88]; % [T(Farenh), % humidity, P(in.Hg)]

[mdata,MeanData,steady_test,phaseSync_test] = loadFiles(source_dir,files_dir,conditions); %load files
[MeanData,StreamData] = loadStreamTripod(mdata,MeanData,source_dir,files_dir); %load streamdata

%Extract variables
collective = mdata.MeanCollective(1);
GR = mdata.GearRatio(1);

%% Process ALL Data
fprintf('\nProcessing data sets...\n')

%Execute remainder of LoadsProcessing.m code
[StreamData,SortedData] = fSortStream(StreamData);

RevData = fRevolutionAvg(SortedData);
AvgData = fTotalAvg(RevData,SortedData,StreamData);


%% Process Steady Data

fprintf('\nProcessing steady-only data...\n')

if exist('steady_test','var')
    [Averages] = runSteady(StreamData,AvgData,steady_test);
end

%% Process Phase Sync
fprintf('\nRunning phase-sync processing...\n')
if exist('phaseSync_test','var')
    offset = input('Offset Angle [deg]: ');
    slew = input('Angle Slew Rate [deg/s]: ');
    [PhaseSync] = runPhaseSync(StreamData,phaseSync_test,offset,Trig,GR,slew);    
end

%% Plotting
close all; clc;

if exist('steady_test','var') && length(steady_test) ~= 1; [f1,f2,f3,f4] = plotSteady(Averages,collective); end
if exist('phaseSync_test','var'); [f5,f6] = plotPhaseSync(PhaseSync,phaseSync_test,loads); end


%% Saving
save_dir = uigetdir();
save(fullfile(save_dir,'Data'),'PhaseSync','Averages');
if exist('f1','var')
    saveas(f1,fullfile(save_dir,f1.Name),'jpg')
    saveas(f2,fullfile(save_dir,f2.Name),'jpg')
    saveas(f3,fullfile(save_dir,f3.Name),'jpg')
    saveas(f4,fullfile(save_dir,f4.Name),'jpg')
end
saveas(f5,fullfile(save_dir,f5.Name),'jpg')
saveas(f6,fullfile(save_dir,f6.Name),'jpg')



