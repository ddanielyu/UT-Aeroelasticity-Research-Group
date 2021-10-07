close all 
figss = dir('*.fig');
figss = {figss.name};
for i = 1 : length(figss)
p = openfig(figss{i});
filename = erase(figss{i},'.fig');
filename = [filename,'.eps'];
saveas(p, filename);
end
close all