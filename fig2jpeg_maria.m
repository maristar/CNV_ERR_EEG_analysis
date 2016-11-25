% this program reads the contents of a folder in figures '.fig' and saves
% them as jpg
files=dir('*Grandaverage*.fig');
for kk=1:length(files); 
    filenames{kk,:}=files(kk,:).name;
end

for m=1:length(files)
    h=open(filenames{m})
    nametemp=filenames{m}(1:end-4)
    saveas(h, nametemp, 'jpeg')
    close(h)
end
