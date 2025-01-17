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

Trig = 10; %rpm spike from nominal rpm to find phase sync trigger

%% Inputs
files_dir = uigetdir(); %directory of .csv data files
motors = input('Single or Dual Motor [s/d]: ','s');
loads = input('Plot with loads? [y/n]: ','s');

%% Change to correct directory for single vs dual motor
if motors == 's'
    cd('Single Motor'); 
    %Check if before or after switch to dual motor
    dual = input('Dual Motor VI? [y/n]: ','s');
elseif motors == 'd'; cd('Dual Motor'); 
elseif isempty(motors); cd('Tripod Pre-Feb 2022'); Trig = 10;
end

source_dir = pwd; %directory of MATLAB scripts
%% Load Data
conditions = [54	29.88]; % [T(Farenh), % humidity, P(in.Hg)]

[mdata,MeanData] = loadFiles(source_dir,files_dir,conditions); %load files
collective = mdata.MeanCollective(1);

if exist('dual','var')
    [MeanData,StreamData] = loadStreamTripod_dual(mdata,MeanData,source_dir,files_dir,dual); %load streamdata
else
    [MeanData,StreamData] = loadStreamTripod(mdata,MeanData,source_dir,files_dir,collective); %load streamdata
end

%Extract variables
GR = mdata.GearRatio(1);

%% Organize Steady and Phase Sync Tests

% Check for phase sync or steady (torque pulse)
fprintf('\nChecking for triggers...\n')
for i = 1:length(StreamData.names)
    nom_IQ = mean(StreamData.IQ{i});
    if contains(StreamData.names{i},'-0') == 0 && exist('dual','var') == 0
        %Phase Sync DID occur
        if exist('phaseSync_test','var'); phaseSync_test{end+1} = StreamData.names{i};
        else; phaseSync_test = {StreamData.names{i}}; end
    
    elseif ismember(1,abs(StreamData.Mz_inner{i} - mean(StreamData.Mz_inner{i})) > Trig) && isempty(motors)
        %Phase Sync DID occur
        if exist('phaseSync_test','var'); phaseSync_test{end+1} = StreamData.names{i};
        else; phaseSync_test = {StreamData.names{i}}; end
        
    elseif ismember(1,abs(StreamData.Mz_inner{i} - mean(StreamData.Mz_inner{i})) > Trig) && exist('dual','var')
        %Phase Sync DID occur
        if exist('phaseSync_test','var'); phaseSync_test{end+1} = StreamData.names{i};
        else; phaseSync_test = mdata.Path; end
            
    else
        %Steady State test
        if exist('steady_test','var'); steady_test{end+1} = StreamData.names{i};
        else; steady_test = {StreamData.names{i}}; end
    end
        
end

%% Process ALL Data
fprintf('\nProcessing data sets...\n')

%Execute remainder of LoadsProcessing.m code
[StreamData,SortedData] = fSortStream(StreamData);

%account for GR
for i = 1:length(StreamData.OMEGA)
    StreamData.OMEGA{i} = StreamData.OMEGA{i}/GR;
    SortedData.cts_outer{i} = SortedData.cts_outer{i}*GR^2;
    SortedData.cps_outer{i} = SortedData.cps_outer{i}*GR^3;
    SortedData.cts_inner{i} = SortedData.cts_inner{i}*GR^2;
    SortedData.cps_inner{i} = SortedData.cps_inner{i}*GR^3;
end

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

if exist('steady_test','var')
    [Averages] = runSteady(StreamData,AvgData,steady_test);
end

%% Process Phase Sync
fprintf('\nRunning phase-sync processing...\n')
if exist('phaseSync_test','var')
    offset = input('Offset Angle [deg]: ');
    if isempty(motors) ~= 1 
        slew = input('Angle Slew Rate [deg/s]: ');
        [PhaseSync] = runPhaseSync(StreamData,phaseSync_test,offset,Trig,GR,slew);
    else 
        [PhaseSync] = runPhaseSync(StreamData,phaseSync_test,offset,Trig,GR);
    end
    
    %clear loads to avoid corrupting plots
    if loads == 'n'
        for i = 1:length(PhaseSync.T_avg)
            PhaseSync.T_avg{i} = NaN;
            PhaseSync.T_err{i} = NaN;
            PhaseSync.Torque_avg{i} = NaN;
            PhaseSync.Torque_err{i} = NaN;
        end
    end
    
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



