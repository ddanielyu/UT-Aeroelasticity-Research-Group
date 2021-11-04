function PL = f_PL(f_c, db)

% CONVERT dB SPECTRA IN 1/3 OCTAVE BANDS INTO PERCIEVED LOUDNESS
% BASED ON HANDBOOK OF NOISE RATINGS

% INPUTS
    % f_c = center frequencies of octave bands, should be 24
    % db = vector of sound in decibels
% OUTPUTS
    % PL = Percieved Loudness, PLdB


load('dBtoSonesTable.mat')

tab_f = sones.freq;
tab_db = sones.db;
tab_s = sones.conversion;
tab_sm = sones.sm;
tab_F = sones.F;

Np = length(f_c); % Number of points
if Np > 24
    f_c = f_c(1:24);
    db = db(1:24);
    N = length(f_c);
    s = ones(N,1);
    
    for i = 1:N
        % Locate closest frequency band to center frequency
        [m,loc_f] = min(abs(f_c(i) - tab_f));
        
        % Get sones value from table
        s(i) = tab_s(round(db(i)) == tab_db, loc_f);
    end
    % Loudness of loudest band
    s_m = max(s);
    
    % Locate F factor
    [m,loc_F] = min(abs(s_m - tab_sm));
    F = tab_F(loc_F);
    
    % Total loudness of sound
    s_t = s_m + F * (sum(s) - s_m);
    
    if s_t>20
        % Percieved loudness
        PL = 32 + 9*log2(s_t); % PLdB
    else
        fprintf('\nCannot calculate percieved loudness, noise too low\n')
        PL = NaN;
    end
    
else
    fprintf('\nCannot calculate percieved loudness, not sufficient number octave bands\n')
    PL = NaN;
end







