function testdata = TestProc_SingleMic(testdate, testletter,micnums, caldata, plots)
% READ AND CONVERT TEST FILES INTO DB VALUES, PLOT EACH MIC
% CMJOHNSON 03/25/2020, EDITED 01/13/2022
% INPUTS
%     caldata                 -> get calibration factors
%     testdate
%     testletter
%     micnums
%     plots = true or false
% OUTPUTS
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



testdata = struct('dbdata',[],'Pdata',[],'Pdata_t',[],'testmag',[],'tvec',[],'wavdata',[],'fs',[],'fvec',[]);
testprefix = [testdate '_test_' testletter ' - 01 Start - '];
% read the test files
for micnum = micnums
    fname = [testprefix num2str(micnum) '.wav'];
    if isfile(fname)
        % load data
        [testdata(micnum).wavdata, testdata(micnum).fs] = audioread(fname);
        x = testdata(micnum).wavdata;
        fs = testdata(micnum).fs;
        N = length(x);
        % process data
        % tvec
        testdata(micnum).tvec = 0: 1/fs: (length(x)-1)/fs;
        
        % fft
        [testdata(micnum).fvec, testdata(micnum).testmag, ~] = ffind_spectrum(fs, x, N, 1);
        mag = testdata(micnum).testmag;
        
        % filter
        window = 'hamming';
        N_avg = 20;
        % 20 avgs with 50% overlap
        [testdata(micnum).fvec_filt, testdata(micnum).testmag_filt, ~] = ffind_spectrum(fs, x, 2*N/N_avg, N_avg, window);
        mag_filt = testdata(micnum).testmag_filt;
        
        % convert to pressure
        % freq domain
        doubling_factor = 1/2;
        testdata(micnum).Pdata = mag * caldata(micnum).calfactor * doubling_factor;
        testdata(micnum).Pdata_filt = mag_filt * caldata(micnum).calfactor * doubling_factor;
        P = testdata(micnum).Pdata;
        P_filt = testdata(micnum).Pdata_filt;
        % time domain
        testdata(micnum).Pdata_t = x * caldata(micnum).calfactor * doubling_factor;
        
        % convert to dB
        Pref = 20E-6; %[Pa]
        testdata(micnum).dbdata = 20*log10(P / Pref);
        testdata(micnum).dbdata_filt = 20*log10(P_filt / Pref);
        
        % OASPL
        testdata(micnum).oaspl = fOverallSPL_time(testdata(micnum).Pdata_t, testdata(micnum).tvec);
        
        fprintf('%s\n',['OASPL = ',num2str(testdata(micnum).oaspl)])
    end
end

if (plots)
    %TIME DOMAIN
    figure(1)
    for micnum = micnums
        subplot(length(micnums),1,micnum)
        plot(testdata(micnum).tvec, testdata(micnum).wavdata)
        axis([0 1 -0.5 0.5]);
        legend(['Mic ' num2str(micnum)]);
    end
    xlabel('Time, s')
    
    
    %FREQ DOMAIN
    figure(11)
    for micnum = micnums
        subplot(length(micnums),1,micnum)
        semilogx(testdata(micnum).fvec, testdata(micnum).dbdata)
        xlim([10^1 10^4]);
        legend(['Mic ' num2str(micnum)]);
    end
    xlabel('Frequency, Hz')
    
    
    figure(13)
    for micnum = micnums
        subplot(length(minnums),1,micnum)
        semilogx(testdata(micnum).fvec, testdata(micnum).dbAdata)
        xlim([10^1 10^4]);
        legend(['Mic ' num2str(micnum)]);
    end
    xlabel('Frequency, Hz')
    
    
end
