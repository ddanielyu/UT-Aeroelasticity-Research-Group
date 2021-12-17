function [MeanData,StreamData] = loadStreamTripod(mdata,MeanData,source_dir,files_dir,collective)
%This function loads all streaming data files from loadFiles.m for the
%2021-2022 Tripod Stand Testing

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
curr3col = 18;       % current 2 column
trigcol = 19;        % analog trigger column
revcol = 20;         % revolution column
rpmcol = 21;

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

switch (rotor)
    case 'Uber'
        MeanData.meancols = ones(length(mdata{:,'MeanCollective'}),1)*collective;
        MeanData.diffcols = mdata{:,'DifferentialCollective'};
        MeanData.zcs = mdata{:,'AxialSpacing'};
        MeanData.phis = mdata{:,'IndexAngle'};
        MeanData.rhos = mdata{:,'rho'};
end

%cd('streaming');   % enter streaming files directory 
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
    StreamData.revolution{k} = [];
    StreamData.trigger{k} = [];
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
    StreamData.curr1{k} = data{:,curr1col};          %W
    StreamData.curr2{k} = data{:,curr2col};          %W
    StreamData.curr3{k} = data{:,curr3col};          %W
    StreamData.revolution{k} = data{:,revcol};       %X
    StreamData.trigger{k} = data{:,trigcol};         %Q           
    StreamData.nrevs{k} = StreamData.revolution{k}(end);
    StreamData.rpm{k} = data{:,rpmcol};
    
    fprintf('%s\n', 'Ok');
    
    StreamData.rho{k} = MeanData.rhos(k);
    
    %% CREATE REV COUNTER
    revnum = 0;
    for i = 1:length(StreamData.encoder{k})-1
        StreamData.revolution{k}(i,1) = revnum;
        if (StreamData.encoder{k}(i) > 359)&& (StreamData.encoder{k}(i+1) < 1)
            revnum = revnum + 1; 
        end
    end
    StreamData.revolution{k}(length(StreamData.encoder{k}),1) = revnum;
    StreamData.nrevs{k} = revnum;
end

cd(source_dir);   % return to original directory

end

