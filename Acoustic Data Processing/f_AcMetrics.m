function [testnames, testdata] = f_AcMetrics(directory)

% READ AND CONVERT TEST FILES INTO 1/3 OCTAVE VALUES AND CALUCALULTE
% METRICS
% CMJOHNSON 10/22/2021

% INPUTS
%     caldata        
%     testdate
%     testletter
%     plots = true or false
%     Pdoubling = true or false

% OUTPUTS
%     testdata
%         .fc = frequency vector (1 x 24)
%         .spl = 1/3 octave values (1 x 24)
%         .oaspl = overall sound pressure level
%         .oaspla = A-weighted sound pressure level
%         .pl = percieved loudness
%         .pnl = percieved noise level

pdir = pwd; 
cd(directory);

%% INPUTS
%Pdoubling
Pdoubling = input('Pressure doubling [y n] ? ', 's');
if (Pdoubling == 'y')
    doubling_factor = 1/2; %PRESSURE DOUBLING AT RIGID SURFACE
else
    doubling_factor = 1;
end

%file names
files = dir('*.wav');
filenames = {files(:).name}';

%choose dates
fprintf('\n%s\n','Test file: ')
dates = unique(extractBefore(filenames,'_'));
fprintf('\n\t%s', 'Loaded test dates are [YYMMDD] : ')
fprintf('%s ',dates{:});
fprintf('\n\t')
testdates = input('Test Date [YYMMDD] : ', 's');
testdates = split(testdates, ' ');

%choose tests
for ii = 1:length(testdates)
    loc = (contains(filenames,'cal'))&(contains(filenames,testdates{ii}));
    caltests = unique(extractBetween(filenames(loc),'test_','_cal'));
    
    loc = (~contains(filenames,'cal'))&(contains(filenames,testdates{ii}));
    tests = unique(extractBetween(filenames(loc),'test_',' -'));
    letters = unique(extractBefore(tests,'_'));
    
    fprintf('\n\t%s', 'Test Date : ')
    fprintf('%s', testdates{ii})
    fprintf('\n\t\t%s', 'Loaded tests are : ')
    fprintf('%s ',letters{:});
    fprintf('\n\t\t')
    testletters{ii} = input('Tests to process : ', 's');
    if strcmp(testletters{ii}, 'all')
        testletters{ii} = letters;
    else
        testletters{ii} = split(testletters{ii}, ' ');
    end

    
    %choose test numbers
    for jj = 1:length(testletters{ii})
        numbers = unique(extractAfter(tests(contains(tests,testletters{ii}{jj})),'_'));
        fprintf('\n\t\t%s', 'Test : ')
        fprintf('%s', testletters{ii}{jj})
        
        fprintf('\n\t\t\t%s', 'Loaded calibration tests are : ')
        fprintf('%s ',caltests{:});
        fprintf('\n\t\t\t')
        calletters{ii}{jj} = input('Calibration test : ', 's');
        calletters{ii}{jj} = split(calletters{ii}{jj}, ' ');
        
        fprintf('\n\t\t\t%s', 'Loaded test numbers are : ')
        fprintf('%s ',numbers{:});
        fprintf('\n\t\t\t')
        testnumbers{ii}{jj} = input('Test numbers to process : ', 's');
        if strcmp(testnumbers{ii}{jj}, 'all')
            testnumbers{ii}{jj} = numbers;
        else 
            testnumbers{ii}{jj} = split(testnumbers{ii}{jj}, ' ');
        end
    end
end

%assemble files names
cnt1=0;
for ii = 1:length(testdates)
    for jj = 1:length(testletters{ii})
        for kk = 1:length(testnumbers{ii}{jj})
            if strcmp(testnumbers{ii}{jj}{kk}, 'all')
                numbers = unique(extractAfter(tests(contains(tests,testletters{ii}{jj})),'_'));
                for ll = 1:length(numbers)    
                    cnt1=cnt1+1;
                    testprefix{cnt1} = [testdates{ii} '_test_' testletters{ii}{jj} '_' numbers{ll} ' - 01 Start - '];
                    calprefix{cnt1} = [testdates{ii} '_test_' calletters{ii}{jj}{1} '_cal - 01 Start - '];
                end
            else
                cnt1=cnt1+1;
                testprefix{cnt1} = [testdates{ii} '_test_' testletters{ii}{jj} '_' testnumbers{ii}{jj}{kk} ' - 01 Start - '];
                calprefix{cnt1} = [testdates{ii} '_test_' calletters{ii}{jj}{1} '_cal - 01 Start - '];
            end
        end
    end
end

%% 
for k = 1:length(testprefix)
% calibration files    
    calname = extractBefore(calprefix{k}, ' -');
    fprintf('\n%s\n',['Calibrating microphones. Calibration file: ',calname])    
    caldata{k} = fCalProc2(calprefix{k});
    
% initialize    
    testdata{k} = struct('name',[],'fc',[],'spl',[],'oaspl',[],'oasplA',[],'oaspl3',[],'pl',[],'pnl',[]);
    testdata{k}.name = extractBefore(testprefix{k},' -');  
    testnames{k} = extractBefore(testprefix{k},' -');  
    fprintf('\n%s\n',['Processing Data. Data file: ', testdata{k}.name])
    for micnum = 1:16
        fname = [testprefix{k} num2str(micnum) '.wav'];
        if isfile(fname)
            fprintf('\t%s',['- Mic ', num2str(micnum),' ... '])
            
% load data            
            [x,fs] = audioread(fname);
            N = length(x);
% process data          
            % tvec
            t = 0: 1/fs: (length(x)-1)/fs;
            
            % fft
            [fvec, mag, ~] = ffind_spectrum(fs, x, N, 1);
            
            % filter
            window = 'hamming';
            N_avg = 20; 
                % 20 avgs with 50% overlap
            [fvec_filt, mag_filt, ~] = ffind_spectrum(fs, x, 2*N/N_avg, N_avg, window);
            
            % convert to pressure
                % freq domain
            P = mag * caldata{k}(micnum).calfactor * doubling_factor; 
            P_filt = mag_filt * caldata{k}(micnum).calfactor * doubling_factor;
                % time domain
            P_t = x * caldata{k}(micnum).calfactor * doubling_factor;
            
                
% octave filtering
                % 1/3
            [testdata{k}(micnum).fc,P3] = fOctaveFilter(fvec,P,3);
            
% convert to dB
            Pref = 20E-6; %[Pa]
            db = 20*log10(P / Pref);
            db_filt = 20*log10(P_filt / Pref);
            db3 = 20*log10(testdata{k}(micnum).ofilt3_Pdata / Pref);
            
% OASPL
            testdata{k}(micnum).oaspl = fOverallSPL_time(P_t, t); 
            testdata{k}(micnum).oasplcheck = fOverallSPL_freq(P);
            testdata{k}(micnum).oasplfilt = fOverallSPL_freq(P_filt);
            testdata{k}(micnum).oaspl3 = fOverallSPL_freq(P3);
            
% A-weighting
            A = fAfilt(fvec);
            dbA = db + A';
            testdata{k}(micnum).oasplA = fOverallSPL_freq(Pref * 10.^(dbA / 20));
            
            A_filt = fAfilt(fvec_filt);
            dbA_filt = db_filt + A_filt';
            testdata{k}(micnum).oasplAfilt = fOverallSPL_freq(Pref * 10.^(dbA_filt / 20));
            
            A3 = fAfilt(testdata{k}(micnum).fc);
            dbA3 = db3 + A3';
            testdata{k}(micnum).oasplA3 = fOverallSPL_freq(Pref * 10.^(dbA3 / 20));
            
% Percieved loudness
            testdata{k}(micnum).PL = f_PL(testdata{k}(micnum).fc, db3);
            
% Percieved noise level
            testdata{k}(micnum).PNL = f_PNL(testdata{k}(micnum).fc, db3);
    
            fprintf('%s\n',['OASPL = ',num2str(testdata{k}(micnum).oaspl)])
        end
    end
    cd(pdir);
    fprintf('\n\t')
end
