%{
This script processes data from Fall 2021 Tripod Stand Testing. Either
steady-state or phase-sync data can be processed.

This processing code is a modification to the "LoadsProcessing.m" code, and
essentially uses only fLoadData

Written By: Matt Asper
Date: 01 Dec 2021

%}

clc; clear; close all;
load('colors.mat');

%% Constants

Torque_Trig = 10; %N-m spike from Nominal Torque to find phase sync trigger

source_dir = pwd; %directory of MATLAB scripts
%% Inputs
files_dir = uigetdir(); %directory of .csv data files
GR = input('Gear Ratio: ');
collective = input('Collective [deg]: ');

%% Load Data
conditions = [54	29.88]; % [T(Farenh), % humidity, P(in.Hg)]

[mdata,MeanData] = loadFiles(source_dir,files_dir,conditions); %load files
[MeanData,StreamData] = loadStreamTripod(mdata,MeanData,source_dir,files_dir,collective); %load streamdata

%% Organize Steady and Phase Sync Tests

% Check for phase sync or steady (torque pulse)
fprintf('\nChecking for triggers...\n')
for i = 1:length(StreamData.names)
    nom_torque = mean(StreamData.Mz_inner{i});
    if ismember(1,abs(StreamData.Mz_inner{i} - nom_torque) > Torque_Trig) 
        %Phase Sync DID occur
        if exist('phaseSync_test','var'); phaseSync_test{end+1} = StreamData.names{i};
        else; phaseSync_test = {StreamData.names{i}}; end
            
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
fprintf('\nProcessing steady-only data...\n')

if exist('steady_test','var')
    [Averages] = runSteady(StreamData,AvgData,steady_test);
end

%% Process Phase Sync
fprintf('\nRunning phase-sync processing...\n')
if exist('phaseSync_test','var')
    offset = input('Offset Angle [deg]: ');
    
    [PhaseSync] = runPhaseSync(StreamData,phaseSync_test,offset,Torque_Trig,GR);
    
end

%% Load Simulink
% plotsim = input('Plot with Sim [y/n]: ','s');
% 
% if plotsim == 'y'
%     fprintf('Select Sim File to Load...\n');
%     [sim_file,sim_path] = uigetfile(strcat('/Users/asper101/Box Sync/For Matt/3rd Year/Electromechanical Modeling/'));
%     load(fullfile(sim_path,sim_file));
%     sim_step = 8;
%     Angle_err.time = Angle_err.time - sim_step;
%     Act_angle.time = Act_angle.time - sim_step;
%     Motor_RPM.time = Motor_RPM.time - sim_step;
%     Q_total.time = Q_total.time - sim_step;
% %     nRev.time = nRev.time - sim_step;
% end
% 
% [~,full_step] = min(abs(nRev.time - sim_step));
% [~,full_step1] = min(abs(nRev.time - sim_step + 0.1));
% nRev_new = nRev.data;

%% Plotting
close all; clc;

if exist('steady_test','var'); plotSteady(Averages,collective); end
if exist('phaseSync_test','var'); plotPhaseSync(PhaseSync); end


%% Saving
% saveas(f1,f1.Name,'jpg')
% saveas(f2,f2.Name,'jpg')



