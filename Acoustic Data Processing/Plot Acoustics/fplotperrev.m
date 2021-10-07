function fplotperrev(RPM,bladenumber)

perrev = RPM/60;
perrevs = perrev*[bladenumber:bladenumber:bladenumber*6]';

YL = get(gca, 'YLim');

plot([perrevs,perrevs],YL,'r-','HandleVisibility','off')

% findpeaks(yd,f, 'MinPeakHeight', 30,'MinPeakDistance',2400/60*2 - 5, 'NPeaks', 30)