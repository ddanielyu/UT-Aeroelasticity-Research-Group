% PROCESS MISC ACOUSTIC DATA 
% CMJOHNSON 06/08/2020
% PROCESS SINGLE TEST ACOUSTIC DATA FILES USING CALIBRATION DATA FILES AND DATA FILES
 
% INPUTS
%     testdate                -> calibration test date / data test date
%     testletter              -> calibration test letter / data test letter
%     calsuffix               -> [] or anything added to end of "cal"
%     plots = true or false   -> plots calibration files
%     filename                -> name of calibration xlsx sheet
%     dirname                 -> location of calibration xlsx sheet
% 
% OUTPUTS
%     caldata
%         .scale              -> max magnitude of wav file at desired calibration frequency
%         .calmag
%         .tvec
%         .wavdata
%         .fs
%         .fvec
%     testdata
%         .fvec
%         .fs
%         .wavdata
%         .tvec
%         .testmag
%         .Pdata [Pa]
%         .Pdata_t [Pa]       -> Pressure in time domain
%         .dbdata
%         .ofilt12_dbdata
%         .ofilt3_dbdata
%         .oaspl

% clear; clc; close all
clear all;
%warning off

%% INPUTS
directory = '\Users\admin-local\Box\Chloe Lab Stuff\2021 Spring Stacked Rotor\Acoustics\zc_15';
directory = '\Users\admin-local\Box\Chloe Lab Stuff\2020 Fall Stacked Rotor\Acoustic Tests';
directory = '../Fall2020AcousticData/';
directory = '/Users/chloe/Box/Chloe Lab Stuff/2020 Fall Stacked Rotor/Acoustic Tests';

addpath(pwd)
inputs.doubling = 'y';

%inputs.date = '201117'; inputs.test = 'f'; inputs.cal = 'a'; 
% inputs.testnum = '5'; %1200 RPM coll:0, i90
% inputs.testnum = '8'; %1200 RPM coll:10, i90
%  inputs.testnum = 'ref'; Background - 0RPM
    
% inputs.date = '201118'; inputs.test = 'b'; inputs.cal = 'a';
%    inputs.testnum = '5'; %1200 RPM coll:0, i90    
%    inputs.testnum = '8'; % 1200 RPM coll:10, i90
%     inputs.testnum = '14 5 6 12 7 11 8 9'; % coll sweep i90 1200 RPM
%inputs.date = '201119'; inputs.test = 'd'; inputs.cal = 'a';
%    inputs.testnum = '5'; %1200 RPM coll:0, i45  
%   inputs.testnum = '8'; % 1200 RPM coll:10, i45

% set environment
env.temp = 20;
env.press= 101325;
env.hr = 30;
far.dist = 72.5;
far.dist_fac = 10;
%% PROCESSING
[testnames, testdata, caldata] = fAcProc(directory,env,far);
fprintf('\n\n%s\n\n', 'Processing done.');

%% PLOT
k=1;
RPM = 1200; %RPM=990;
bladenumber = 2;
micnum = 3;

figure(1); hold on;
plot(testdata{k}(micnum).tvec-7.4,testdata{k}(micnum).Pdata_t)


figure(2); % vs elevtion
plot(linspace(-90,90,16),[testdata{k}.oaspl])

figure(3);  % spectrum
%semilogx(tetestdata{k}(micnum).fvec, testdata{k}(micnum).dbdata); hold on
%semilogx(testdata{k}(micnum).ofilt12_fvec, testdata{k}(micnum).ofilt12_dbdata,'linewidth',1.1)
semilogx(testdata{k}(micnum).fvec_filt, testdata{k}(micnum).dbdata_filt,'k', 'displayname','Experiment'); hold on;
semilogx(testdata{k}(micnum).fvec_filt, testdata{k}(micnum).dBbb, 'displayname','Experiment filt'); hold on;
semilogx(testdata{k}(micnum).fvec_filt, testdata{k}(micnum).dBtl, 'displayname','Experiment filt'); hold on;

%semilogx(testdata{k}(micnum).fvec_filt, testdata{k}(micnum).dBfarA,'k', 'displayname','Experiment far'); hold on;

xlim([10^1 2*10^4]);
ylim([0 80])
xlabel('Frequency, Hz')
ylabel('SPL, dB \Delta f = 1Hz')
% fplotperrev(RPM,bladenumber)

%% save
% exp = table(); %data structure for saving
% % make vector for specific mic from multiple tests
% for i = 1:length(testdata)
% exp.names{i} = [testdata{i}(1).name];
% exp.db(i,1) = [testdata{i}(micnum).oaspl];
% exp.dbA(i,1) = [testdata{i}(micnum).oasplA];
% exp.dBbb(i,1) = [testdata{i}(micnum).oasplbb];
% exp.dBtl(i,1) = [testdata{i}(micnum).oaspltl];
% exp.dbfar(i,1) = [testdata{i}(micnum).oasplfar];
% exp.dbfarbb(i,1) = [testdata{i}(micnum).oasplfarbb];
% exp.dbfartl(i,1) = [testdata{i}(micnum).oasplfartl];
% exp.dbfarA(i,1) = [testdata{i}(micnum).oasplfarA];
% exp.dbfarAbb(i,1) = [testdata{i}(micnum).oasplfarAbb];
% exp.dbfarAtl(i,1) = [testdata{i}(micnum).oasplfarAtl];
% end
% 
% 
% figure(4);
% plot(exp.db); hold on;
% plot(exp.dBbb)
% plot(exp.dBtl)
addpath('\Users\admin-local\Box\Chloe Lab Stuff\2021 Spring Stacked Rotor\Results')
% load('zc15_ac_v2.mat')

% for i = length(data)+1:length(data)+1+length(testdata)
for i = 1:length(testdata)
%     j = i -length(data);
j = i;
    data{i}(1).name = [testdata{j}(1).name];
for micnum = 1:14
data{i}(micnum).f = [testdata{j}(micnum).fvec_filt];
data{i}(micnum).db = [testdata{j}(micnum).dbdata_filt];
data{i}(micnum).dbA = [testdata{j}(micnum).dbAdata_filt];
data{i}(micnum).oaspl = [testdata{j}(micnum).oaspl];
data{i}(micnum).oasplA = [testdata{j}(micnum).oasplA];
data{i}(micnum).oasplAbb = [testdata{j}(micnum).oasplAbb];
data{i}(micnum).oasplAtl = [testdata{j}(micnum).oasplAtl];
data{i}(micnum).dBbb = [testdata{j}(micnum).oasplbb];
data{i}(micnum).dBtl = [testdata{j}(micnum).oaspltl];
data{i}(micnum).dbfar = [testdata{j}(micnum).oasplfar];
data{i}(micnum).dbfarbb = [testdata{j}(micnum).oasplfarbb];
data{i}(micnum).dbfartl = [testdata{j}(micnum).oasplfartl];
data{i}(micnum).dbfarA = [testdata{j}(micnum).oasplfarA];
data{i}(micnum).dbfarAbb = [testdata{j}(micnum).oasplfarAbb];
data{i}(micnum).dbfarAtl = [testdata{j}(micnum).oasplfarAtl];
% data{i}(micnum).testbb = [testdata{j}(micnum).oaspl_tonal];
% data{i}(micnum).testtl = [testdata{j}(micnum).oaspl_bb];
end
end