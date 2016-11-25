% Make statistical analysis 07-12-2011 on grandaverages
% takes the resultscor and it calculates the intra-hemispheric
% connectivity. revised 16.12.2011
clear all 
close all

%% GO grandaverage
cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program')
DOI='D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program';

nchan=input('Number of channels ');
numRight=[1,2,20:29]; % Cental electrodes excluded
numLeft=[4, 6:15,17];
N={'FR2', 'FZ2', 'FCZ', 'CZ1', 'FZA', 'FZ1', 'FL1', 'FL3', 'FL5','CL3', 'CL1', 'CL5', 'PL5', 'PL1', 'PL3','PZC', 'O1', 'PZP', 'OZ', 'O2', 'PR4', 'PR2', 'PR6', 'CR2', 'CZ2', 'CR6', 'CR4', 'FR6', 'FR4'};

triggerlist={'go', 'nogo'};
textmeasuresall={'cor','pcor','DTFdelta','DTFtheta','DTFalpha','DTFbeta','DTFgamma'};
for qq=1:length(textmeasuresall)
    textmeasure=textmeasuresall{qq} % what is the measure
    %% for trigger GO
    trigger=triggerlist{1};
    [av, XYZ] = grandaverager_m(DOI, trigger, textmeasure, nchan)
    cd(DOI)
    mkdir(textmeasure)
    cd(textmeasure)
    resultsALL.(textmeasure).(trigger).grandav=av;
    clear grandaverage_go
    %save resultsALL resultsALL
    %% for trigger NOGO take the grandaverage
    trigger=triggerlist{2};
    grandaverage_nogo=zeros(nchan);
    [grandaverage_nogo, XYZ]=grandaverager_m(DOI, trigger, textmeasure, nchan);
    cd(DOI)
    cd(textmeasure)
    resultsALL.(textmeasure).(trigger).grandav=grandaverage_nogo;
    clear grandaverage_nogo grandaverage_go
    cd(DOI)
    save resultsALL resultsALL
    %% Left - Right Assymetry
    s=N;
    resultsALL.XYZ=XYZ;
    corgofR=zeros(nchan); %initialize 
    corgofL=zeros(nchan); %initialize
    corgoall=zeros(nchan);
    %% for go and nogo separately
    for kkk=1:length(triggerlist)
       trigger=triggerlist{kkk};
       mcor=resultsALL.(textmeasure).(trigger).grandav;
       %nchan=length(XYZ);
       corgofR(numRight, numRight)= mcor(numRight,numRight);
       corgofL(numLeft, numLeft)= mcor(numLeft,numLeft); 
       corgoall=[corgofR+corgofL]./2;
       thr(kkk)=(1/2)*(max(max(squeeze(max(corgoall)))))% to miso tou megistou
       grand_lines(textmeasure,corgofR,corgofL, XYZ, s, DOI, max(thr), numLeft, numRight,trigger);
       resultsALL.(textmeasure).(trigger).conn_L=corgofL;
       resultsALL.(textmeasure).(trigger).conn_R=corgofR;
       resultsALL.(textmeasure).(trigger).conn_thr=thr;
       save resultsALL resultsALL -v7.3
       clear corgof corgofL thr corgoall
       corgof=zeros(nchan); %initialize 
       corgofL=zeros(nchan); %initialize
       corgoall=zeros(nchan);
   end
end

cd(DOI)
save resultsALL resultsALL