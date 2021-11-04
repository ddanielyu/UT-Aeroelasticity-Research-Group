function PNL = f_PNL(f_c, db)

% CONVERT dB SPECTRA IN 1/3 OCTAVE BANDS INTO PERCIEVED NOISE LEVEL
% BASED ON HANDBOOK OF NOISE RATINGS

% INPUTS
    % f_c = center frequencies of octave bands, should be 24
    % db = vector of sound in decibels
% OUTPUTS
    % PNL = Percieved Noise Level, PNdB

load('PNLConversionTable.mat');

Np = length(f_c); % Number of points
if Np > 24
    f_c = f_c(1:24);
    db = db(1:24);
    N = length(f_c);
    noys = ones(N,1);
    
    for i = 1:N
        % Sound pressure level in the band
        L = db(i);
        
        % Locate octave band
        [m,loc_f] = min(abs(f_c(i) - pnl_conv.f));
        L1 = pnl_conv.L1(loc_f);
        L2 = pnl_conv.L2(loc_f);
        L3 = pnl_conv.L3(loc_f);
        L4 = pnl_conv.L4(loc_f);
        Lc = pnl_conv.Lc(loc_f);
        
        if L1<=L && L<L2
            M1 = pnl_conv.M1(loc_f);
            noys(i) = 0.1 * 10^(M1 * (L-L1));
            
        elseif L2<=L && L<L3
            M2 = pnl_conv.M2(loc_f);
            noys(i) = 10^(M2 * (L-L3));
            
        elseif L3<=L && L<Lc
            M3 = pnl_conv.M3(loc_f);
            noys(i) = 10^(M3 * (L-L3));
            
        elseif Lc<=L && L<=150
            M4 = pnl_conv.M4(loc_f);
            noys(i) = 10^(M4 * (L-L4));
        end
    end
    % Noys in the band with greatest noys
    n_m = max(noys);
 
    % Total noy value
    n_t = n_m + 0.15 * (sum(noys) - n_m);
    

    % Percieved noise level
    PNL = 40 + 33.22*log10(n_t); % PNdB

else
    fprintf('\nCannot calculate percieved noise level, not sufficient number octave bands\n')
    PNL = NaN;
end
