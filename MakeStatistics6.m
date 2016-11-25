% Make statistical analysis 07-12-2011 on grandaverages
% takes the resultscor and extracts the measure then calculates the
% grandaveage among many channels. Then it calculates the intra-hemispheric
% connectivity. revised 16.12.2011, used for savoury 20-12-2011
clear all 
close all

%% GO grandaverage
cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\Second component')
DOI='D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\Second component';
% cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program')
% DOI='D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program';
cd GO
files=dir('*');
for kk=1:length(files); 
    filenames{kk,:}=files(kk,:).name;
end
ND=kk; clear kk
nchan=input('Number of channels ')
textmeasuresall={'cor','pcor','DTFdelta','DTFtheta','DTFalpha','DTFbeta','DTFgamma'};
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
        cd(stemp)
        load resultscor
        resultsALL.XYZ=resultscor.XYZ;
        %resultsALL.s=resultscor.s;
        switch textmeasure
            case 'pcor', meantemp=squeeze(mean(resultscor.(textmeasure).result2,3));
            case 'cor', meantemp=squeeze(mean(resultscor.(textmeasure).result1,3));
            case {'DTFdelta','DTFtheta','DTFalpha','DTFbeta'}, meantemp=resultscor.resultsDTF.(textmeasure).meangamma; %for %dtf
            case 'DTFgamma1', textmeasurex='DTFgamma', meantemp=resultscor.resultsDTF.(textmeasurex).meangamma;
        end
        grandaverage=meantemp+grandaverage;
        currentNtrial=resultscor.num_epochs;
        Ntrial=Ntrial+currentNtrial;
          chan_strength_norm=resultscor.(textmeasure).chan_strength_norm;
        
          switch textmeasure
            case {'pcor','cor', 'DTFdelta','DTFtheta','DTFalpha','DTFbeta'}, chan_appearance_norm=resultscor.(textmeasure).chan_appearance; %
            case 'DTFgamma1', textmeasurex='DTFgamma', chan_appearance_norm=resultscor.(textmeasure).chan_appearance_norm; %meantemp=resultscor.resultsDTF.(textmeasurex).meangamma;
        end
        
          %_norm; %% mporei na einai chan_appearance_norm
        [br_areas]=brainar(chan_strength_norm, chan_appearance_norm);
        resultscor.(textmeasure).br_areas=br_areas;
        allcouples(kkm,:)=resultscor.(textmeasure).couples.couple_conn_values;
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
    files2=dir('*');
    for kk=1:length(files); 
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
            case {'DTFdelta','DTFtheta','DTFalpha','DTFbeta'}, meantemp=resultscor.resultsDTF.(textmeasure).meangamma; %for %dtf
            case 'DTFgamma1', textmeasurex='DTFgamma', meantemp=resultscor.resultsDTF.(textmeasurex).meangamma;
        end
        grandaverage_nogo=meantemp+grandaverage_nogo;
       currentNtrial=resultscor.num_epochs;
        Ntrial=Ntrial+currentNtrial;
        chan_strength_norm=resultscor.(textmeasure).chan_strength_norm;
        
        switch textmeasure
            case {'pcor','cor', 'DTFdelta','DTFtheta','DTFalpha','DTFbeta'}, chan_appearance_norm=resultscor.(textmeasure).chan_appearance; %
            case 'DTFgamma1', textmeasurex='DTFgamma', chan_appearance_norm=resultscor.(textmeasure).chan_appearance_norm; %meantemp=resultscor.resultsDTF.(textmeasurex).meangamma;
        end
        
        
        
        [br_areas]=brainar(chan_strength_norm, chan_appearance_norm);
        resultscor.(textmeasure).br_areas=br_areas;
        allcouples_nogo(kkm,:)=resultscor.(textmeasure).couples.couple_conn_values;
        clear meantemp resultscor
    cd ..   
    end % end files for Nogo %%%%%%%%%%%%%%%%%%%%
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
numRight=[1,2,20:29]; % Cental electrodes excluded
numLeft=[4, 6:15,17]; % ok!
% s=resultsALL.s;

nchan=length(XYZ);
corgof=zeros(nchan);
corgofL=zeros(nchan);
corgoall=zeros(nchan);
corgof(numRight, numRight)=mcor(numRight,numRight);
corgofL(numLeft, numLeft)=mcor(numLeft,numLeft); 
corgoall=(corgof+corgofL)./2;
thr=0.5*(max(squeeze(max(corgoall))));% to miso tou megistou
%titi=1-0.05^(1/(Ntrial_go-1)); % as Haliday 1995
grand_lines(textmeasure,corgof,corgofL, XYZ, s, DOI, thr, numLeft, numRight,trigger);


%% for NOGO 
trigger=triggerlist{2};
s=N;
XYZ=resultsALL.XYZ;
%textmeasure='DTFtheta';
mcor=resultsALL.(textmeasure).nogo;

numRight=[1,2,20:29]; % Cental electrodes excluded
numLeft=[4, 6:15,17];
% s=resultsALL.s;

nchan=length(XYZ);
corgof=zeros(nchan);
corgofL=zeros(nchan);
corgoall=zeros(nchan);
corgof(numRight, numRight)=mcor(numRight,numRight);
corgofL(numLeft, numLeft)=mcor(numLeft,numLeft); 
corgoall=(corgof+corgofL)./2;
thr=0.8*(max(squeeze(max(corgoall))));% to 80% tou megistou
titi=1-0.05^(1/(Ntrial_nogo-1)); % as Haliday 1995
grand_lines(textmeasure,corgof,corgofL, XYZ, s, DOI, thr, numLeft, numRight,trigger);



end


    for jj=3:18, [p(jj), h(jj)]=ttest2(measure1(jj,:), measure1_nogo(jj,:)); end  %% file-wise statistics
    figure; plot(h, '*');
    for kk=1:nchan, [pc(kk), hc(kk)]=ttest2(measure1(3:end,kk), measure1_nogo(3:end,kk)); end %% channel-wise statistics
    figure; plot(hc, '*'); 
cd(DOI)
save resultsALL resultsALL