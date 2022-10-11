function [MeanData,StreamData] = loadStreamTripod(mdata,MeanData,source_dir,files_dir,collective,SF,GR,rpm_smoothing)
%This function loads all streaming data files from loadFiles.m for the
%2021-2022 Tripod Stand Testing
% FEB 14, 2022 UPDATE: accommodates tripod 220212 Labview (streamed phase 3
% and IQ, no analog trig, no rpm col)

cd(files_dir);

%% DECLARE QUANTITIES IN STREAMING DATA FILE

Fxocol = 1;          % Fx_outer column
Fyocol = 2;          % Fy_outer column
Fzocol = 3;          % Fz_outer column
Mxocol = 4;          % Mx_outer column
Myocol = 5;          % My_outer column
Mzocol = 6;          % Mz_outer column
Fxicol = 7;          % Fx_inner column
Fyicol = 8;          % Fy_inner column
Fzicol = 9;          % Fz_inner column
Mxicol = 10;         % Mx_inner column
Myicol = 11;         % My_inner column
Mzicol = 12;         % Mz_inner column
axcol = 13;          % mag Ax column
aycol = 14;          % mag Ay column    
enccol = 15;         % encoder angle column
curr1col = 16;       % current 1 column
curr2col = 17;       % current 2 column
buscol = 18;         % bus current column
curr3col = 19;       % current 3 column
IQcol   = 20;        % IQ column
revcol = 21;         % revolution column


%% FIND NAMES OF STREAMING FILES AND ASSIGN OPERATING VARIABLES

StreamData.names = mdata{:,'Path'};

fprintf('\n%s\n', 'Checking streaming files ...');

% remove rows corresponding to files that dont exist
row2rm = false(length(StreamData.names),1);   % vector of row numbers to remove
for ii = 1:length(StreamData.names)
    if ~isfile(StreamData.names{ii})
        row2rm(ii) = true;
        fprintf('\t%s\n', ['Missing file : ' StreamData.names{ii}]);
    end
end
mdata(row2rm,:) = [];
StreamData.names(row2rm,:) = [];

% find nominal RPM : actual RPM within closest 5 RPM. 898 -> 900, 901 -> 900
MeanData.RPMs = (round(mdata{:,'RPM'}/5))*5;  


%% LOAD AND PROCESS STREAMING FILES
rotor = 'Uber';

% switch (rotor)
%     case 'Uber'
%         MeanData.meancols = ones(length(mdata{:,'MeanCollective'}),1)*collective;
%         MeanData.rhos = mdata{:,'rho'};
% end

% cd('Streaming');   % enter streaming files directory 
[nfiles, ~] = size(mdata);

fprintf('\n%s\n', 'Reading streaming files');

for k = 1:nfiles
    StreamData.Fx_outer{k} = [];
    StreamData.Fy_outer{k} = [];
    StreamData.Fz_outer{k} = [];
    StreamData.Mx_outer{k} = [];
    StreamData.My_outer{k} = [];
    StreamData.Mz_outer{k} = [];
    StreamData.Fx_inner{k} = [];
    StreamData.Fy_inner{k} = [];
    StreamData.Fz_inner{k} = [];
    StreamData.Mx_inner{k} = [];
    StreamData.My_inner{k} = [];
    StreamData.Mz_inner{k} = [];
    StreamData.ax{k} = [];
    StreamData.ay{k} = [];
    StreamData.encoder{k} = [];
    StreamData.curr1{k} = [];
    StreamData.curr2{k} = [];
    StreamData.curr3{k} = [];
    StreamData.bus_curr{k} = [];
    StreamData.revolution{k} = [];
    StreamData.nrevs{k} = [];
    
    fprintf('\t%s', ['- ' StreamData.names{k} ' ... ']);

    data = readtable(StreamData.names{k});
    StreamData.Fx_outer{k} = data{:,Fxocol};         %A
    StreamData.Fy_outer{k} = data{:,Fyocol};         %B
    StreamData.Fz_outer{k} = data{:,Fzocol};    %C
    StreamData.Mx_outer{k} = data{:,Mxocol};         %D
    StreamData.My_outer{k} = data{:,Myocol};         %E
    StreamData.Mz_outer{k} = data{:,Mzocol};      %F
    StreamData.Fx_inner{k} = data{:,Fxicol};         %G
    StreamData.Fy_inner{k} = data{:,Fyicol};         %H
    StreamData.Fz_inner{k} = data{:,Fzicol};         %I
    StreamData.Mx_inner{k} = data{:,Mxicol};         %J
    StreamData.My_inner{k} = data{:,Myicol};         %K
    StreamData.Mz_inner{k} = data{:,Mzicol};         %L
    StreamData.ax{k} = data{:,axcol};                %M
    StreamData.ay{k} = data{:,aycol};                %N
    StreamData.encoder{k} = abs(data{:,enccol});          %W
    StreamData.unwrap{k} = unwrap_az(StreamData.encoder{k}); %unwrap the azimuth
    StreamData.curr1{k} = data{:,curr1col};          %W
    StreamData.curr2{k} = data{:,curr2col};          %W
    StreamData.curr3{k} = -1*(StreamData.curr1{k} + StreamData.curr2{k});          %W
    StreamData.IQ{k} = parkClarke(StreamData.curr1{k},StreamData.curr2{k},StreamData.curr3{k}); 
    StreamData.bus_curr{k} = data{:,buscol};          %W
    StreamData.rpm{k} = (StreamData.unwrap{k}(2:end) - StreamData.unwrap{k}(1:end-1))*SF/360*60;
    StreamData.rpm{k} = fcleanup(StreamData.rpm{k}, 'smoothdata', 'loess', rpm_smoothing);
%     StreamData.revolution{k} = data{:,revcol} - data{1,revcol};       %X
%     StreamData.nrevs{k} = StreamData.revolution{k}(end);
    
    fprintf('%s\n', 'Ok');
    
    StreamData.rho{k} = MeanData.rhos(k);
    
    %% CREATE REV COUNTER
    revnum = 0;
    for i = 1:length(StreamData.encoder{k})-1
        StreamData.revolution{k}(i,1) = revnum;
%         if (StreamData.encoder{k}(i) > 359)&& (StreamData.encoder{k}(i+1) < 1)
%             revnum = revnum + 1; 
%         end
        if StreamData.encoder{k}(i) >= StreamData.encoder{k}(i+1)
            revnum = revnum + 1; 
        end
    end
    StreamData.revolution{k}(length(StreamData.encoder{k}),1) = revnum;
    StreamData.nrevs{k} = revnum;
end

cd(source_dir);   % return to original directory

end

%% Functions

function [unwrapped_az] = unwrap_az(wrapped_az)
%This function unwraps the encoder azimuth readings from 0-360 deg to 0-inf
%
%Written By: Matt Asper
%Date: 09 August 2022

unwrapped_az(1) = 0;

for i = 2:length(wrapped_az)
    
    if wrapped_az(i) >= wrapped_az(i-1)
        unwrapped_az(i) = unwrapped_az(i-1) + wrapped_az(i) - wrapped_az(i-1);
    elseif wrapped_az(i) < wrapped_az(i-1)
        unwrapped_az(i) = unwrapped_az(i-1) + wrapped_az(i) + (360 - wrapped_az(i-1));
    end
        
end

end

function Q = parkClarke(A,B,C)
%This function calculates the quadrature (Q) current from A,
%B, and C phase currents by assuming all current produced from controller
%is IQ
%{
Inputs: 
A,B,C -- phase currents (column vectors)

Outputs:
Q -- Quadrature current (column vector)
%}

for i = 1:length(A)
    Q(i) = sqrt(2/3)*norm([A(i) B(i) C(i)]);
end
Q = Q';

end


function [filt_Vs] = fcleanup(Vs, method, arg, winlen, fl, fs, fr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to clean up time signal using either a filter or a smoothing
% algorithm
%
% INPUTS:
%        Vs = signal time series
%        method = filtering/ smoothing algorithm
%        arg = number of harmonics or number of data points for smoothing method
%        winlen = number of elements at the ends to zero out
%        fl = low-pass cut-off frequency, Hz
%        fs = sampling frequency, Hz
%        fr = 1/rev frequency, Hz
%
% OUTPUTS:
%        filt_Vs = filtered signal time series
%
% sirohi 191129
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 7) % First low pass the signal and then band stop 'arg' harmonics
    
    % Define a windowing function to make the ends of the signal = 0
    trad  = linspace(0,1,winlen)*pi;
    winfn = [0.5-cos(trad)*0.5 ones(1,length(Vs)-2*winlen) cos(trad)*0.5+0.5]';
    wVs   = Vs .*winfn;
    
    % Lowpass filter
    filt_Vs = lowpass(wVs, fl, fs , 'StopbandAttenuation', 80, 'Steepness', 0.95);
    
    % Bandstop filter
    for ii = 1:arg
        
        filt_Vs = bandstop(filt_Vs, [ii*fr-5 ii*fr+5], fs);  % Remove ii/rev +/- 5Hz
    
    end
    
else
    
    if strcmp(method, 'movmean') % Moving average
       
        eval(['filt_Vs = ' method '(Vs, ' num2str(arg) ');']);

    elseif strcmp(method, 'smoothdata') % Smoothing filter defined in matlab
        
        eval(['filt_Vs = ' method '(Vs, arg , ' num2str(winlen) ');']);
  
    end
    
end

% End of function

end

