% Import and filter OFC patients, 21-1-2012, MLS
clear all 
close all

cd('D:\OFC');
DOI='D:\OFC';
files=dir('*.raw');
for kk=1:length(files)/2; 
    filenames{kk,:}=files(kk+(length(files)/2),:).name;
end
ND=kk;
clear kk
thismoment=datestr(now);
for jj=1:length(thismoment); if ( thismoment(jj)==':' || thismoment(jj)==' '); thismoment(jj)='-';end; end
clear jj files
tic
for kkm=1:ND
    disp(filenames{kkm}(1:end-5));
    disp(kkm)
        cd([DOI '\SETS\']);
        %stemp=([filenames{kkm}(1:end-5) '-' thismoment])
        stemp=([filenames{kkm}(1:end-5)]);
        mkdir(stemp)
        cd(stemp)
        dirtemp=pwd;
        %% (1) Load the set with the raw dataset
        cd(DOI);
        EEG=pop_readegi_m(filenames{kkm});
        pop_saveset(EEG,'filename', filenames{kkm}(1:end-5),'filepath',dirtemp, 'savemode', 'resave');
        Fs=EEG.srate;
        DCoffset_removal_21_10_2011_b_2012;
        EEG = pop_chanedit(EEG, 'load', {'D:\RIKSHOSPITALET\CNV_RIKS\RAW DATASETS\EGI129-HydroCel-Original.spf', 'filetype', 'sfp'});
        stemp2=[filenames{kkm}(1:end-5) '_filt'];
        EEG = pop_saveset(EEG, 'filename', stemp2 , 'check', 'on',  ...
        'filepath',dirtemp); % no pop-up D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\'
        clear EEG dirtemp stemp2 stemp data_filt data legen1 k h d 
        
        %pop_saveset(EEG,'filename', filenames{kkm}(1:end-5),'filepath',dirtemp, 'savemode', 'resave')
        % start filter
end


