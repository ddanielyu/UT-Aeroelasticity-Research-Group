% ACOUSTIC METRICS PROCESSING
% CHLOE JOHNSON 10/22/2021

clear; clc; close all
warning off

%% INPUTS
directory = '\Users\cmj2855\Box\Chloe Lab Stuff\2021 Spring Stacked Rotor\Acoustics\zc_15';

%% PROCESSING
[testnames, testdata] = f_AcMetrics(directory);
fprintf('\n\n%s\n\n', 'Processing done.');
