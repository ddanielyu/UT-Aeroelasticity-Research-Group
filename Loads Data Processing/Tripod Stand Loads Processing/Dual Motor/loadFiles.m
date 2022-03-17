function [mdata,MeanData,steady_test,phaseSync_test] = loadFiles(source_dir,files_dir,conditions)
%This scripts loads all the data files for Tripod Stand 2021-2022 data
%processing

%% Inputs
T_F = input('Temperature [F]: ');


cd(files_dir);
Files = dir('*.csv');
FileName = {Files.name};
dates = unique(extractBefore(FileName,'_'));

fprintf('\n\t%s', 'Loaded test dates are [YYMMDD] : ')
fprintf('%s ',dates{:});
fprintf('\n\t')
testdates = input('[YYMMDD] : ', 's');
testdates = split(testdates, ' ');

jj = 1;

Files = FileName(contains(FileName,'mean') & contains(FileName,testdates{jj}));
letters = extractBetween(Files,'test_',' mean');
fprintf('\n\t%s', 'Test Date : ')
fprintf('%s', testdates{jj})
fprintf('\n\t%s', 'Loaded test letters are : ')
fprintf('%s ',letters{:});
fprintf('\n\t')

steady_letters{jj} = input('Steady test letters : ','s'); 
steady_letters{jj} = split(steady_letters{jj},' ');


phase_sync_letters{jj} = input('Phase sync letters : ','s');
phase_sync_letters{jj} = split(phase_sync_letters{jj},' ');
testletters{jj} = [steady_letters{jj};phase_sync_letters{jj}];

%% Processing
cnt = 0;
for ii = 1:length(testletters{jj})
    if strcmp(testletters{jj}{ii}, 'all')
        Files = FileName(contains(FileName,'mean') & contains(FileName,testdates{jj}));
        letters = extractBetween(Files,'test_',' mean');
        for kk = 1:length(letters)
            cnt = cnt+1;
            testnames{cnt} = [testdates{jj}, '_test_' letters{kk}, ' mean'];
        end        
    else
        testletters{jj}{ii} = ['_test_' testletters{jj}{ii} ' mean'];
        cnt = cnt+1;
        testnames{cnt} = [testdates{jj}, testletters{jj}{ii}];
    end
end


TF = contains(FileName,testnames);
MeanData.names = FileName(TF);
mdata = table(); % assemble a table with all the mean data
for im = 1:length(MeanData.names) % need to fix reading multiple mean files
    fprintf('\n%s\t\n\t', ['Loading mean data file : ' MeanData.names{im} ' ...']);

    opts = detectImportOptions(MeanData.names{im});
    opts.DataLines = [2 Inf];
    opts.VariableNamesLine = 1;
    if im>1
        opts.VariableNames = vnames;
    end
    if im == 1
        MeanData.data{im} = readtable(MeanData.names{im}, opts, 'ReadVariableNames', true);
        vnames = MeanData.data{im}.Properties.VariableNames;
    else
        MeanData.data{im} = readtable(MeanData.names{im}, opts, 'ReadVariableNames', false);
    end
    
    T = (T_F - 32)*5/9 + 273.15; % [Kelvin]
    P = conditions(2)*101325/29.9212; % [Pa]
    R_air = 287.05; % INDIVIDUAL GAS CONSTANT
    rho = P/R_air/T;
    MeanData.data{im}.rho = rho*ones(height(MeanData.data{im}),1);
    
    ExVar = startsWith(MeanData.data{im}.Properties.VariableNames,'ExtraVar');
    if sum(ExVar) > 0
        MeanData.data{im}.ExtraVar1 = [];
    end
        
    mdata = vertcat(mdata, MeanData.data{im});
end    

%% Organize Steady and Phase Sync Tests
steady_test = [];
phaseSync_test = [];

ij = 1; %steady test counter
ji = 1; %phase sync test counter

for i = 1:length(mdata{:,'Path'})
    %steady
    for ii = 1:length(steady_letters)
        if isempty(steady_letters{1}); break; end
        
        if contains(mdata{i,'Path'},strcat('test_',steady_letters{ii}))
            steady_test{ij} = mdata{i,'Path'};
            ij = ij + 1;
        end
    end
    
    %phase-sync
    for jj = 1:length(phase_sync_letters)
        if isempty(phase_sync_letters{1}); break; end
        
        if contains(mdata{i,'Path'},strcat('test_',phase_sync_letters{jj}))
            phaseSync_test{ji} = mdata{i,'Path'};
            ji = ji + 1;
        end
    end
end

% for i = 1:length(phaseSync_tests)
%     
% end

cd(source_dir);

end

