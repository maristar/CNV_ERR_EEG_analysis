 % Make statistical analysis 07-12-2011 on grandaverages, 27-2-2012
% takes the resultscor and extracts the measure then calculates the
% grandaveage among many channels. Then it calculates the intra-hemispheric
% connectivity. revised 16.12.2011, used for savoury 20-12-2011
clear all 
close all

%% GO grandaverage
% cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program')
% DOI='D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program';

%DOI='D:\OFC\ANALYZED_DATASETS_cube\All_interval_OFC\'
%DOI='D:\OFC\ANALYZED_DATASETS_cube\Second_component_OFC\'
%DOI='D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS_cube\NEW\all_interval_controls\'
% DOI='D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS_cube\NEW\Second_component_controls\'
DOI='/Users/mstavrin/Documents/MATLAB/CNV/ANALYZED_DATASETS/'
cd(DOI)
cd GO
files=dir('*ICA*');
for kk=1:length(files); 
    filenames{kk,:}=files(kk,:).name;
end
ND=kk; clear kk
nchan=input('Number of channels ')
textmeasuresall={'cor','pcor','DTFdelta','DTFtheta','DTFalpha','DTFbeta','DTFgamma1'};
for qq=1:length(textmeasuresall)
    textmeasure=textmeasuresall{qq};
    grandaverage=zeros(nchan);
    grandaverage_nogo=zeros(nchan);
    Ntrial=0;
    for kkm=3:ND
        currentNtrial=0;
        cd(DOI)
        cd GO
        %disp(filenames(kkm));
        disp(kkm)
        stemp=([filenames{kkm}]);
        disp(stemp)
        cd(stemp)
        load resultscor
        resultsALL.XYZ=resultscor.XYZ;
        %resultsALL.s=resultscor.s;
        switch textmeasure
            case 'pcor', meantemp=squeeze(mean(resultscor.(textmeasure).result2,3));
            case 'cor', meantemp=squeeze(mean(resultscor.(textmeasure).result1,3));
            case {'DTFdelta','DTFtheta','DTFalpha','DTFbeta','DTFgamma1'}, meantemp=resultscor.resultsDTF.(textmeasure).meangamma; %for %dtf
        end
        grandaverage=meantemp+grandaverage;
        currentNtrial=resultscor.num_epochs;
        Ntrial=Ntrial+currentNtrial;
        clear meantemp resultscor
    end
    grandaverage=grandaverage./(ND-2);
    cd(DOI)
    mkdir(textmeasure)
    cd(textmeasure)
    resultsALL.(textmeasure).go=grandaverage;
    clear kkm 
    %% NO GO
    cd(DOI)
    cd NOGO
    files2=dir('*ICA*');
    for kk=1:length(files2); 
        filenames2{kk,:}=files2(kk,:).name;
    end
    ND=kk; clear kk; Ntrial=0;
    for kkm=3:ND %% start files for NOGO
        cd(DOI)
        cd NOGO
        disp(kkm)
        currentNtrial=0;
        stemp=([filenames2{kkm}]);
        cd(stemp)
        load resultscor
        switch textmeasure
            case 'pcor', meantemp=squeeze(mean(resultscor.(textmeasure).result2,3));
            case 'cor', meantemp=squeeze(mean(resultscor.(textmeasure).result1,3));
            case {'DTFdelta','DTFtheta','DTFalpha','DTFbeta','DTFgamma1'}, meantemp=resultscor.resultsDTF.(textmeasure).meangamma; %for %dtf
        end
        grandaverage_nogo=meantemp+grandaverage_nogo;
       currentNtrial=resultscor.num_epochs;
        Ntrial=Ntrial+currentNtrial;
        clear meantemp resultscor
    cd ..   
    end % end files for Nogo
    grandaverage_nogo=grandaverage_nogo./(ND-2);
    Ntrial_nogo=Ntrial./(ND-2);
    cd(DOI)
    cd(textmeasure)
    resultsALL.(textmeasure).nogo=grandaverage_nogo;
    
    %measure1_nogo=squeeze(measure1_nogo);
    XYZ=resultsALL.XYZ;  
    N={'FR2', 'FZ2', 'FCZ', 'CZ1', 'FZA', 'FZ1', 'FL1', 'FL3', 'FL5','CL3', 'CL1', 'CL5', 'PL5', 'PL1', 'PL3','PZC', 'O1', 'PZP', 'OZ', 'O2', 'PR4', 'PR2', 'PR6', 'CR2', 'CZ2', 'CR6', 'CR4', 'FR6', 'FR4'};
    resultsALL.(textmeasure).s=N;
    clear resultscor

triggerlist={'go', 'nogo'};
trigger=triggerlist{1}
s=N;
XYZ=resultsALL.XYZ;
mcor=resultsALL.(textmeasure).go;
numRight=[1,21,22,23,24,26,27,28,29]; % 22-2-2012\\\ old: numRight=[1,2,20:29]; numLeft=[4, 6:15,17]; % ok!% Cental electrodes excluded
numLeft=[7:15]; % ok!
% s=resultsALL.s;
name='Grandaverage';
nchan=length(XYZ);
corgof_go=zeros(nchan);
corgofL_go=zeros(nchan);
corgoall_go=zeros(nchan);
corgof_go(numRight, numRight)=mcor(numRight,numRight);
corgofL_go(numLeft, numLeft)=mcor(numLeft,numLeft); 
corgoall_go=(corgof_go+corgofL_go)./2;
thr_go=0.5*(max(squeeze(max(corgoall_go))));% to miso tou megistou
%titi=1-0.05^(1/(Ntrial_go-1)); % as Haliday 1995



%% for NOGO 
trigger=triggerlist{2};
s=N;
XYZ=resultsALL.XYZ;
%textmeasure='DTFtheta';
mcor=resultsALL.(textmeasure).nogo;

numRight=[1,21,22,23,24,26,27,28,29]; % 22-2-2012\\\ old: numRight=[1,2,20:29]; numLeft=[4, 6:15,17]; % ok!% Cental electrodes excluded
numLeft=[7:15];
% s=resultsALL.s;

nchan=length(XYZ);
corgof_nogo=zeros(nchan);
corgofL_nogo=zeros(nchan);
corgoall_nogo=zeros(nchan);
corgof_nogo(numRight, numRight)=mcor(numRight,numRight);
corgofL_nogo(numLeft, numLeft)=mcor(numLeft,numLeft); 
corgoall_nogo=(corgof_nogo+corgofL_nogo)./2;
thr_nogo=0.8*(max(squeeze(max(corgoall_nogo))));% to 80% tou megistou
titi=1-0.05^(1/(Ntrial_nogo-1)); % as Haliday 1995

thr=min(thr_go,thr_nogo);

grand_lines(textmeasure,name,corgof_go,corgofL_go, XYZ, s, DOI, thr, numLeft, numRight,triggerlist{1});
grand_lines(textmeasure,name, corgof_nogo,corgofL_nogo, XYZ, s, DOI, thr, numLeft, numRight,triggerlist{2});

resultsALL.(textmeasure).corgof_go=corgof_go;
resultsALL.(textmeasure).corgof_nogo=corgof_nogo;
resultsALL.(textmeasure).corgofL_nogo=corgofL_nogo;
resultsALL.(textmeasure).corgofL_go=corgofL_go;
resultsALL.(textmeasure).corgoall_nogo=corgoall_nogo;
resultsALL.(textmeasure).corgoall_go=corgoall_go;
resultsALL.(textmeasure).thr_nogo=thr_nogo;
resultsALL.(textmeasure).thr_go=thr_go;
resultsALL.(textmeasure).thr=thr;

end

cd(DOI)
save resultsALL resultsALL

% for qq=1:length(textmeasuresall)
%     c

%clear corgof_go corgofL_go corgofL_nogo corgof_nogo thr_nogo thr_go