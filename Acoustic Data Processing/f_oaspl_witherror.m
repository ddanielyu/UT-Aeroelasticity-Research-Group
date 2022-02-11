function [dbVs,dbX, dbx] = f_oaspl_witherror(fs, x, NFFT, Nav, win, params)
%--------------------------------------------------------------------------
% function to return mean and standard deviation of oaspl of signal
% Vs must be uniformly sampled, i.e., dt = constant
% divides the input signal into blocks and averages the spectra
%
% INPUTS :
% fs   - sampling frequency, Hz
% x    - signal row vector
% NFFT - number of elements to use for FFT (must be even)
% Nav  - number of averages (must be odd)
% win  - window type - omit if no window required
%
% OUTPUTS :
% avg_oaspl
% sd_oaspl
% johnson 211025 modified from ffind_spectrum
%--------------------------------------------------------------------------
%%
% pad the end of the vectors if length of the vector is odd
if ( mod(length(x),2) == 1 )
    x = [x x(end)];
end

N = length(x);
n = 0:N-1;          % time indices
tvec =  n/fs;       % time vector
k = 0:NFFT-1;       % frequency indices
f = k *fs/NFFT;     % frequency vector

if nargin==6
    eval(['winfn = ' win '(NFFT);']);
    winfn = winfn .*sqrt(NFFT/sum(winfn.^2)); % normalize window for P= 1 W
else
    winfn = rectwin(NFFT);
end

%% create matrix with NFFT rows and Nav columns, each containing a block of x
xblock = zeros(NFFT, Nav);
tblock = zeros(NFFT, Nav);

if Nav>1
    Noffset = floor((N-NFFT) /(Nav-1));
else
    Noffset = 0;
end

if Noffset==0
    Nav = 1;
end

Noverlap = (1 -Noffset/NFFT)*100;
% disp(['Block size: ' num2str(NFFT) ' Averages: ' num2str(Nav) ' Overlap: '...
%     num2str(Noverlap) '%']);

for iav = 1:Nav
    idx1 = (iav-1) *Noffset +1;
    idx2 = idx1 +NFFT -1;
    xblock(:,iav) = x(idx1:idx2);
    tblock(:,iav) = tvec(idx1:idx2);
    xblock(:,iav) = xblock(:,iav).*winfn;
end

%% find the DFTs, scale, average
DFTXblock2 = fft(xblock);              % matrix containing the DFT of the columns of xblock
DFTXblock = DFTXblock2(1:NFFT/2, :);   % single-sided DFT
magXblock = abs(DFTXblock);             % absolute value
magVs2 = 2 /NFFT *magXblock;           % magnitude at each frequency
powerVs2 = 2* magXblock.^2/ NFFT^2;      % power in each bin
sd_magVs = std(magVs2,0,2);
sd_mag = mean(sd_magVs);              % in .wav units
magX = mean(magVs2,2);                % average of all the blocks
powerX = mean(powerVs2,2);

DFTx = fft(x); 
DFTx = DFTx(1:N/2);
magx = abs(DFTx);
magx = 2 /N *magx;                   % magnitude at each frequency
k = 0:N-1;       % frequency indices
fx = k *fs/N;     % frequency vector

fvec = f(1:NFFT/2);                    % generate frequency vector
fvecx = fx(1:N/2);

% figure()
% for i =1:Nav
% subplot(Nav,1,i)
% plot(tblock(:,i),xblock(:,i))
% xlim([0,10])
% end

%% process data
%convert to pressure
% freq domain
PVs = magVs2 .* params.cal .* params.doub;
PX = magX .* params.cal .* params.doub;
Px = magx .* params.cal .* params.doub;

% convert to dB
Pref = 20E-6; %[Pa]
dbVs = 20*log10(PVs ./ Pref);
dbX = 20*log10(PX ./ Pref);
dbx = 20*log10(Px ./ Pref);

% OASPL
oasplVs = fOverallSPL_freq(PVs);
oasplX = fOverallSPL_freq(PX);
oasplx = fOverallSPL_freq(Px);

% A-weighting
% A = fAfilt(fvec);
% dbAVs = dbVs + A';
% dbAX = dbX + A';
% oasplAVs = fOverallSPL_freq(Pref * 10.^(dbAVs / 20));
% oasplAX = fOverallSPL_freq(Pref * 10.^(dbAX / 20));


end