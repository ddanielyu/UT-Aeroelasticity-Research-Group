function [testnames, testdata, caldata] = fAcProc(directory,env,far,varargin)
% READ AND CONVERT TEST FILES INTO DB VALUES, PLOT EACH MIC
% CMJOHNSON 05/15/2020
% INPUTS
%     env: .temp, .press, .rh - temperature (C), pressure (Pa), rel. humidity(%)
%     caldata                 -> get calibration factors
%     testdate
%     testletter
%     plots = true or false
%     Pdoubling = true or false
%     varargin = structure of inputs [ doubling, date,test,cal,testnum]
% OUTPUTS
%     testdata
%         .fvec               -> frequency vector (1 x 240000)
%         .fs                 -> sampling frequency (48000 Hz)
%         .wavdata            -> data in .wav file (480000 x 1)
%         .tvec               -> time vector (1 x 240000)
%         .testmag            -> magnitudes of .wav file in freq. domain
%                                (240000 x 1)
%         .Pdata [Pa]         -> pressure magnitudes in freq.
%                                domain (240000 x 1)
%         .dbdata             -> pressure magnitudes in freq. domain
%                                converted to dB (240000 x 1)
%         .Pdata_t [Pa]       -> pressure in time domain (240000 x 1)
%
%         .ofilt12_fvec
%         .ofilt12_Pdata
%         .ofilt12_dbdata
%
%         .ofilt3_fvec
%         .ofilt3_Pdata
%         .ofilt3_dbdata
%
%         .oaspl
%         .dbdataA
%         .oasplA
pdir = pwd;
cd(directory);
if length(varargin) >0
    getin = false;
    inputs = varargin{1};
else
    getin = true;
end
%% INPUTS

if getin
    Pdoubling = input('Pressure doubling [y n] ? ', 's');
else
    Pdoubling = inputs.doubling;
end

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
if getin
    testdates = input('Test Date [YYMMDD] : ', 's');
else
    testdates = inputs.date;
end
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
    if getin
        testletters{ii} = input('Tests to process : ', 's');
    else
        testletters{ii} = inputs.test;
    end
    testletters{ii} = split(testletters{ii}, ' ');
    
    %choose test numbers
    for jj = 1:length(testletters{ii})
        numbers = unique(extractAfter(tests(contains(tests,testletters{ii}{jj})),'_'));
        fprintf('\n\t\t%s', 'Test : ')
        fprintf('%s', testletters{ii}{jj})
        
        fprintf('\n\t\t\t%s', 'Loaded calibration tests are : ')
        fprintf('%s ',caltests{:});
        fprintf('\n\t\t\t')
        if getin
            calletters{ii}{jj} = input('Calibration test : ', 's');
        else
            calletters{ii}{jj} = inputs.cal;
        end
        calletters{ii}{jj} = split(calletters{ii}{jj}, ' ');
        
        fprintf('\n\t\t\t%s', 'Loaded test numbers are : ')
        fprintf('%s ',numbers{:});
        fprintf('\n\t\t\t')
        if getin
            testnumbers{ii}{jj} = input('Test numbers to process : ', 's');
        else
            testnumbers{ii}{jj} = inputs.testnum;
        end
        
        if strcmp(testnumbers{ii}{jj},'all')
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
    testdata{k} = struct('name',[],'oaspl',[],'oasplA',[],'dbAdata',[],'dbdata',[],'Pdata',[],'Pdata_t',[],'testmag',[],'tvec',[],'wavdata',[],'fs',[],'fvec',[],'ofilt12_fvec',[],'ofilt12_Pdata',[],'ofilt12_dbdata',[],'ofilt3_fvec',[],'ofilt3_Pdata',[],'ofilt3_dbdata',[]);
    testdata{k}.name = extractBefore(testprefix{k},' -');  
    testnames{k} = extractBefore(testprefix{k},' -');  
    fprintf('\n%s\n',['Processing Data. Data file: ', testdata{k}.name])
%     fprintf('%s\n',['      SPL    SPLtl   SPLbb   SPLfar    SPLfartl  SPLfarbb  SPLfarA    SPLfarAtl  SPLfarAbb'])
    for micnum = 1:16
        fname = [testprefix{k} num2str(micnum) '.wav'];
        if isfile(fname)
%             fprintf('\t%s',['- Mic ', num2str(micnum),' ... '])
            
% load data            
            [x, testdata{k}(micnum).fs] = audioread(fname);
%             x = testdata{k}(micnum).wavdata;
            fs = testdata{k}(micnum).fs;
            N = length(x);
% process data          
            % tvec
            t = 0: 1/fs: (length(x)-1)/fs;
%             t = testdata{k}(micnum).tvec;
            
            % fft
            [f, mag, ~] = ffind_spectrum(fs, x, N, 1);
%             f = testdata{k}(micnum).fvec;
%             mag = testdata{k}(micnum).testmag;
            
            % filter
            window = 'hamming';
            N_avg = 20; % N_avg = 20; % 20 avgs with 50% overlap
            [testdata{k}(micnum).fvec_filt, testdata{k}(micnum).testmag_filt, ~] = ffind_spectrum(fs, x, 2*N/N_avg, N_avg, window);
            mag_filt = testdata{k}(micnum).testmag_filt;
            
            % convert to pressure
                % freq domain
            P = mag * caldata{k}(micnum).calfactor * doubling_factor; 
            testdata{k}(micnum).Pdata_filt = mag_filt * caldata{k}(micnum).calfactor * doubling_factor;
%             P = testdata{k}(micnum).Pdata;
            P_filt = testdata{k}(micnum).Pdata_filt;            
                % time domain
            testdata{k}(micnum).Pdata_t = x * caldata{k}(micnum).calfactor * doubling_factor;
            P_t = testdata{k}(micnum).Pdata_t;
            
            % phase averaging 
%             encname = [testprefix{k} num2str(24) '.wav'];
%             enc = audioread(encname);
%             if max(enc) > 1e-3
%                 [testdata{k}(micnum).P_sort, P_tonal, P_bb, fs_new] = fPhaseAvg(P_t, enc);
%                 N = length(P_tonal);
%                 % freq domain
%                 [testdata{k}(micnum).fvec_tonal, P_spectra_tonal, ~] = ffind_spectrum(fs_new, P_tonal, N,1);
%                 [testdata{k}(micnum).fvec_bb, P_spectra_bb, ~] = ffind_spectrum(fs_new,P_bb, 2*N/N_avg, N_avg, window);
%                  % db
%                 Pref = 20E-6; %[Pa]
%                 testdata{k}(micnum).db_tonal = 20*log10(P_spectra_tonal / Pref);
%                 testdata{k}(micnum).db_bb = 20*log10(P_spectra_bb / Pref);
%                 % oaspl
%                 testdata{k}(micnum).oaspl_tonal = fOverallSPL_freq(P_spectra_tonal);
%                 testdata{k}(micnum).oaspl_bb = fOverallSPL_freq(P_spectra_bb);   
%             end
            
            
            % octave filtering
                % 1/12
            %[testdata{k}(micnum).ofilt12_fvec,testdata{k}(micnum).ofilt12_Pdata] = fOctaveFilter(f,P,12);
                % 1/3
            %[testdata{k}(micnum).ofilt3_fvec,testdata{k}(micnum).ofilt3_Pdata] = fOctaveFilter(f,P,3);
            
            % convert to dB
            Pref = 20E-6; %[Pa]
%             testdata{k}(micnum).dbdata = 20*log10(P / Pref);
            testdata{k}(micnum).dbdata_filt = 20*log10(P_filt / Pref);
            %testdata{k}(micnum).ofilt12_dbdata = 20*log10(testdata{k}(micnum).ofilt12_Pdata / Pref);
            %testdata{k}(micnum).ofilt3_dbdata = 20*log10(testdata{k}(micnum).ofilt3_Pdata / Pref);
            
            % filter spectrum
            % remove_peaks(bpf in Hz, df in Hz, npeaks, pressure^2, inner%, outer%)
            [pmsbb,pmstl] = remove_peaks(40,1,35,P_filt.^2,15,25);
            testdata{k}(micnum).dBbb = 10*log10(pmsbb/Pref^2);
            testdata{k}(micnum).dBtl = 10*log10(pmstl/Pref^2);
            
            % propagation
            testdata{k}(micnum).dBfar = propagate(testdata{k}(micnum).fvec_filt',testdata{k}(micnum).dbdata_filt,far.dist,far.dist_fac,env.temp,env.press,env.hr);
            testdata{k}(micnum).dBfarbb = propagate(testdata{k}(micnum).fvec_filt',testdata{k}(micnum).dBbb,far.dist,far.dist_fac,env.temp,env.press,env.hr);
            testdata{k}(micnum).dBfartl = propagate(testdata{k}(micnum).fvec_filt',testdata{k}(micnum).dBtl,far.dist,far.dist_fac,env.temp,env.press,env.hr);
            
            
            % A-weighting
%             A = fAfilt(f);
%             testdata{k}(micnum).dbAdata = testdata{k}(micnum).dbdata + A';
            A_filt = fAfilt(testdata{k}(micnum).fvec_filt);
            testdata{k}(micnum).dbAdata_filt = testdata{k}(micnum).dbdata_filt + A_filt';
            testdata{k}(micnum).dBAbb = testdata{k}(micnum).dBbb+ A_filt';
            testdata{k}(micnum).dBAtl = testdata{k}(micnum).dBtl+ A_filt';
            testdata{k}(micnum).dBfarA = testdata{k}(micnum).dBfar+ A_filt';
            testdata{k}(micnum).dBfarAbb = testdata{k}(micnum).dBfarbb+ A_filt';
            testdata{k}(micnum).dBfarAtl = testdata{k}(micnum).dBfartl+ A_filt';
            
            % OASPL
            testdata{k}(micnum).oaspl = fOverallSPL_time(testdata{k}(micnum).Pdata_t, t);
            testdata{k}(micnum).oasplA = fOverallSPL_freq(Pref * 10.^(testdata{k}(micnum).dbAdata_filt / 20));
            testdata{k}(micnum).oasplbb = fOverallSPL_freq(sqrt(pmsbb));
            testdata{k}(micnum).oaspltl = fOverallSPL_freq(sqrt(pmstl));
            testdata{k}(micnum).oasplAbb = fOverallSPL_freq(Pref * 10.^(testdata{k}(micnum).dBAbb / 20));
            testdata{k}(micnum).oasplAtl = fOverallSPL_freq(Pref * 10.^(testdata{k}(micnum).dBAtl / 20));
            testdata{k}(micnum).oasplfar = fOverallSPL_freq(Pref * 10.^(testdata{k}(micnum).dBfar / 20));
            testdata{k}(micnum).oasplfarbb = fOverallSPL_freq(Pref * 10.^(testdata{k}(micnum).dBfarbb / 20));
            testdata{k}(micnum).oasplfartl = fOverallSPL_freq(Pref * 10.^(testdata{k}(micnum).dBfartl / 20));
            testdata{k}(micnum).oasplfarA = fOverallSPL_freq(Pref * 10.^(testdata{k}(micnum).dBfarA / 20));
            testdata{k}(micnum).oasplfarAbb = fOverallSPL_freq(Pref * 10.^(testdata{k}(micnum).dBfarAbb / 20));
            testdata{k}(micnum).oasplfarAtl = fOverallSPL_freq(Pref * 10.^(testdata{k}(micnum).dBfarAtl / 20));
    
            
            
%             fprintf('%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f \n',testdata{k}(micnum).oaspl,testdata{k}(micnum).oaspltl,testdata{k}(micnum).oasplbb,testdata{k}(micnum).oasplfar,testdata{k}(micnum).oasplfartl, testdata{k}(micnum).oasplfarbb,testdata{k}(micnum).oasplfarA,testdata{k}(micnum).oasplfarAtl, testdata{k}(micnum).oasplfarAbb)
        end
    end
    %cd(pdir);
    fprintf('\n\t')
%% VISUALIZE    
%     worv = input('Visualize (v) test data ? ', 's');
worv = 'n';
    switch worv
        case 'v'
            
            %FREQ DOMAIN
            figure(11)
            for micnum = 1:8
                subplot(8,1,micnum)
                semilogx(testdata{k}(micnum).fvec, testdata{k}(micnum).dbdata,testdata{k}(micnum).fvec, testdata{k}(micnum).dbAdata)
                xlim([10^1 10^4]);
                legend(['Mic ' num2str(micnum)]);
            end
            xlabel('Frequency, Hz')
            sgtitle([testdate '_test_' testletter])
            
            figure(12)
            for micnum = 9:16
                subplot(8,1,micnum-8)
                semilogx(testdata{k}(micnum).fvec, testdata{k}(micnum).dbdata,testdata{k}(micnum).fvec, testdata{k}(micnum).dbAdata)
                xlim([10^1 10^4]);
                legend(['Mic ' num2str(micnum)]);
            end
            xlabel('Frequency, Hz')
            sgtitle([testdate ' test ' testletter])
            
        otherwise
    end
end
end