function [MeanData,StreamData] = loadStreamTripod(mdata,MeanData,source_dir,files_dir)
%This function loads all streaming data files from loadFiles.m for the
%2021-2022 Tripod Stand Dual Motor Testing
% March 15, 2022 UPDATE: accommodates dual motor tripod 220315 Labview 

cd(files_dir);
% nCurrents = input('Number of current sensors [3/5]: ');

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
enc1col = 15;         % ROTOR servo encoder angle column
enc2col = 16;         % ROTOR follower encoder angle column
indexcol = 24;         % ROTOR index angle column
curr1col = 17;       % motor 1 current 1 column
curr2col = 18;       % motor 1 current 2 column
buscol = 19;         % bus current column
curr4col = 20;       % motor 2 current 1 column
curr5col = 21;       % motor 2 current 2 column
revcol1 = 22;         % ROTOR servo revolution column
revcol2 = 23;         % ROTOR follower revolution column

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
MeanData.ServoRPMs = (round(mdata{:,'ServoRPM'}/5))*5;  
MeanData.FollowerRPMs = (round(mdata{:,'FollowerRPM'}/5))*5;  


%% LOAD AND PROCESS STREAMING FILES
rotor = 'Uber';

switch (rotor)
    case 'Uber'
        MeanData.meancols = mdata{:,'MeanCollective'};
        MeanData.uppercols = mdata{:,'UpperCollective'};
        MeanData.lowercols = mdata{:,'LowerCollective'};
        MeanData.zcs = mdata{:,'AxialSpacing'};
        MeanData.phis = mdata{:,'IndexAngle'};
        MeanData.rhos = mdata{:,'rho'};
end

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
    StreamData.encoder1{k} = [];
    StreamData.encoder2{k} = [];
    StreamData.index{k} = [];
    StreamData.curr1{k} = [];
    StreamData.curr2{k} = [];
    StreamData.curr3{k} = [];
    StreamData.curr4{k} = [];
    StreamData.curr5{k} = [];
    StreamData.curr6{k} = [];
    StreamData.bus_curr{k} = [];
    StreamData.revolution1{k} = [];
    StreamData.revolution2{k} = [];
    StreamData.nrev1{k} = [];
    StreamData.nrev2{k} = [];
    
    fprintf('\t%s', ['- ' StreamData.names{k} ' ... ']);

    data = readtable(StreamData.names{k});
    StreamData.Fx_outer{k} = data{:,Fxocol};         %A
    StreamData.Fy_outer{k} = data{:,Fyocol};         %B
    StreamData.Fz_outer{k} = data{:,Fzocol} * -1;    %C
    StreamData.Mx_outer{k} = data{:,Mxocol};         %D
    StreamData.My_outer{k} = data{:,Myocol};         %E
    StreamData.Mz_outer{k} = data{:,Mzocol} * -1;      %F
    StreamData.Fx_inner{k} = data{:,Fxicol};         %G
    StreamData.Fy_inner{k} = data{:,Fyicol};         %H
    StreamData.Fz_inner{k} = data{:,Fzicol};         %I
    StreamData.Mx_inner{k} = data{:,Mxicol};         %J
    StreamData.My_inner{k} = data{:,Myicol};         %K
    StreamData.Mz_inner{k} = data{:,Mzicol};         %L
    StreamData.ax{k} = data{:,axcol};                %M
    StreamData.ay{k} = data{:,aycol};                %N
    StreamData.encoder1{k} = data{:,enc1col};                %N
    StreamData.encoder2{k} = data{:,enc2col};                %N
    StreamData.index{k} = data{:,indexcol};                %N
    StreamData.curr1{k} = data{:,curr1col};          %W
    StreamData.curr2{k} = data{:,curr2col};          %W
    StreamData.curr3{k} = -1*(StreamData.curr1{k} + StreamData.curr2{k});          %W
    StreamData.curr4{k} = data{:,curr4col};
    StreamData.curr5{k} = data{:,curr5col};
    StreamData.curr6{k} = -1*(StreamData.curr4{k} + StreamData.curr5{k});
    StreamData.IQ1{k} = parkClarke(StreamData.curr1{k},StreamData.curr2{k},StreamData.curr3{k}); 
    StreamData.IQ2{k} = parkClarke(StreamData.curr4{k},StreamData.curr5{k},StreamData.curr6{k});

    
    StreamData.bus_curr{k} = data{:,buscol};          %W
    StreamData.revolution1{k} = data{:,revcol1};       %X
    StreamData.revolution2{k} = data{:,revcol2};       %X
    
    %bias revs incase they start at a non-zero number
    StreamData.revolution1{k} = StreamData.revolution1{k} - StreamData.revolution1{k}(1);
    StreamData.revolution2{k} = StreamData.revolution2{k} - StreamData.revolution2{k}(1);
    
    StreamData.nrev1{k} = StreamData.revolution1{k}(end);
    StreamData.nrev2{k} = StreamData.revolution2{k}(end);
    
    %account for mismatch between revs
    if StreamData.nrev1{k} > StreamData.nrev2{k}
        StreamData.nrev1{k} = StreamData.nrev2{k}; 
    elseif StreamData.nrev2{k} > StreamData.nrev1{k}
        StreamData.nrev2{k} = StreamData.nrev1{k};
    end
    
    fprintf('%s\n', 'Ok');
    
    StreamData.rho{k} = MeanData.rhos(k);
    
    %"zero" encoders so each encoder starts at 0 deg
    StreamData.encoder1{k} = StreamData.encoder1{k} - StreamData.encoder1{k}(1);
    StreamData.encoder2{k} = StreamData.encoder2{k} - StreamData.encoder2{k}(1);
    
end

cd(source_dir);   % return to original directory

end

