% Make statistical analysis 07-12-2011, rev 23-12-2011
clear all 
close all
% cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\Second component')
% DOI='D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\Second component';

%doi='D:\OFC\SETS_correct_locations\'; %'D:\RIKSHOSPITALET\CNV_RIKS\RAW DATASETS\After_eye_art_removal'
DOI='D:\OFC\ANALYZED_DATASETS\'; 
cd(DOI)

thismoment=datestr(now);
for jj=1:length(thismoment); if ( thismoment(jj)==':' || thismoment(jj)==' '); thismoment(jj)='-';end; end

cd GO %% for the GO triggers
files=dir('*');
for kk=1:length(files); 
    filenames{kk,:}=files(kk,:).name; % to filenames einai 3 parapanw.
end
ND=kk;
nchan=input('Number of channels')
textmeasuresall={'cor','pcor','DTFdelta','DTFtheta','DTFalpha','DTFbeta','DTFgamma1'}
%% PART 1. MAKE THE STATSFINAL WITH ALL VARIABLES INSIDE 
clear jj kk files 
for qq=1:length(textmeasuresall)
    textmeasure=textmeasuresall{qq}
    allcouples_nogo=0;
    allcouples=0;
    temp=0;
    grandaverage=zeros(nchan);
    grandaverage_nogo=zeros(nchan);
    for kkm=1:(ND-2)
        cd(DOI)
        cd GO
        %disp(filenames(kkm));
        disp(kkm)
        stemp=([filenames{kkm+2}]);
        cd(stemp)
        load resultscor
        % measure 0, grandaverages & results 
        switch textmeasure
            case 'pcor', meantemp=squeeze(mean(resultscor.(textmeasure).result2,3));
            case 'cor', meantemp=squeeze(mean(resultscor.(textmeasure).result1,3));
            case {'DTFdelta','DTFtheta','DTFalpha','DTFbeta'}, meantemp=resultscor.resultsDTF.(textmeasure).meangamma; %for %dtf
            case 'DTFgamma1', meantemp=resultscor.resultsDTF.(textmeasure).meangamma;
        end
        grandaverage=meantemp+grandaverage;
        measure0(kkm,:,:)=meantemp;
        % measure .1. chan strength
        chan_strength_norm=resultscor.(textmeasure).chan_strength_norm;
        % chan appearance switch 
        switch textmeasure
            case {'pcor','cor','DTFtheta','DTFalpha','DTFbeta'}, chan_appearance_norm=resultscor.(textmeasure).chan_appearance; %
            case {'DTFgamma1', 'DTFdelta'}, chan_appearance_norm=resultscor.(textmeasure).chan_appearance_norm; %meantemp=resultscor.resultsDTF.(textmeasurex).meangamma;
        end
        
%         [br_areas]=brainar(chan_strength_norm, chan_appearance_norm);
%         resultscor.(textmeasure).br_areas=br_areas; save resultscor resultscor -v7.3
        measure1(kkm, :,:)=resultscor.(textmeasure).chan_strength_norm;
        
        % measure .2. couple strength
        temp=resultscor.(textmeasure).couples.couple_conn_values;
        allcouples=temp+allcouples;
        measure2(kkm,:,:)=temp;
        
        % measure .3. areas strength
        measure3_areas(kkm,1,:)=resultscor.(textmeasure).br_areas.FR;
        measure3_areas(kkm,2,:)=resultscor.(textmeasure).br_areas.FL;
        measure3_areas(kkm,3,:)=resultscor.(textmeasure).br_areas.FZ;
        measure3_areas(kkm,4,:)=resultscor.(textmeasure).br_areas.CR;
        measure3_areas(kkm,5,:)=resultscor.(textmeasure).br_areas.CL;
        measure3_areas(kkm,6,:)=resultscor.(textmeasure).br_areas.CZ;
        measure3_areas(kkm,7,:)=resultscor.(textmeasure).br_areas.PR;
        measure3_areas(kkm,8,:)=resultscor.(textmeasure).br_areas.PL;
        measure3_areas(kkm,9,:)=resultscor.(textmeasure).br_areas.PZ;
        measure3_areas(kkm,10,:)=resultscor.(textmeasure).br_areas.OZ;
        clear meantemp resultscor temp
    end
    allcouples=allcouples./(ND-2);
    grandaverage=grandaverage./(ND-2);
    %resultsstatsALL.(textmeasure).go=grandaverage;
    measure2=squeeze(measure2);
    measure1=squeeze(measure1);
    measure3_areas=squeeze(measure3_areas);
    clear chan_appearance_norm chan_strength_norm temp stemp meantemp
    
    %% NO GO
    cd(DOI)
    cd NOGO
    files2=dir('*');
    for kk=1:length(files2); 
        filenames2{kk,:}=files2(kk,:).name;
    end
    ND=kk;
    
    for kkm=1:(ND-2)
        cd(DOI)
        cd NOGO
        %disp(filenames2(kkm));
        disp(kkm)
        stemp=([filenames2{kkm+2}]);
        cd(stemp)
        load resultscor
        %measure 0 results.
        switch textmeasure
            case 'pcor', meantemp=squeeze(mean(resultscor.(textmeasure).result2,3));
            case 'cor', meantemp=squeeze(mean(resultscor.(textmeasure).result1,3));
            case {'DTFdelta','DTFtheta','DTFalpha','DTFbeta'}, meantemp=resultscor.resultsDTF.(textmeasure).meangamma; %for %dtf
            case 'DTFgamma1', meantemp=resultscor.resultsDTF.(textmeasure).meangamma;
        end
        grandaverage_nogo=meantemp+grandaverage_nogo;
        measure0_nogo(kkm,:,:)=meantemp;
        
        % measure .1. chan strength
        chan_strength_norm=resultscor.(textmeasure).chan_strength_norm;
        
         % chan appearance switch 
        switch textmeasure
            case {'pcor','cor', 'DTFtheta','DTFalpha','DTFbeta'}, chan_appearance_norm=resultscor.(textmeasure).chan_appearance; %
            case {'DTFgamma1', 'DTFdelta'}, textmeasurex='DTFgamma', chan_appearance_norm=resultscor.(textmeasure).chan_appearance_norm; %meantemp=resultscor.resultsDTF.(textmeasurex).meangamma;
        end
        
       
%         [br_areas]=brainar(chan_strength_norm, chan_appearance_norm);
%         resultscor.(textmeasure).br_areas=br_areas;
%         save resultscor resultscor -v7.3
        measure1_nogo(kkm, :,:)=resultscor.(textmeasure).chan_strength_norm;
        
        % measure .2. couple strength
        temp=resultscor.(textmeasure).couples.couple_conn_values;
        allcouples_nogo=temp+allcouples_nogo;
        measure2_nogo(kkm,:,:)=temp;
        
        % measure .3. areas strength
        measure3_areas_nogo(kkm,1,:)=resultscor.(textmeasure).br_areas.FL;
        measure3_areas_nogo(kkm,2,:)=resultscor.(textmeasure).br_areas.FR;
        measure3_areas_nogo(kkm,3,:)=resultscor.(textmeasure).br_areas.CZ;
        measure3_areas_nogo(kkm,4,:)=resultscor.(textmeasure).br_areas.CL;
        measure3_areas_nogo(kkm,5,:)=resultscor.(textmeasure).br_areas.CR;
        measure3_areas_nogo(kkm,6,:)=resultscor.(textmeasure).br_areas.PZ;
        measure3_areas_nogo(kkm,7,:)=resultscor.(textmeasure).br_areas.PL;
        measure3_areas_nogo(kkm,8,:)=resultscor.(textmeasure).br_areas.PR;
        measure3_areas_nogo(kkm,9,:)=resultscor.(textmeasure).br_areas.OZ;
        measure3_areas_nogo(kkm,10,:)=resultscor.(textmeasure).br_areas.FZ;
        clear meantemp temp % to resultscor to kratame gia na paroume ta XYZ
        
       
               
    end
    cd(DOI)
    allcouples_nogo=allcouples_nogo./(ND-2);
    grandaverage_nogo=grandaverage_nogo./(ND-2);
    %resultsstatsALL.(textmeasure).nogo=grandaverage_nogo;
    measure1_nogo=squeeze(measure1_nogo);
    measure3_areas_nogo=squeeze(measure3_areas_nogo);
    measure2_nogo=squeeze(measure2_nogo);
    XYZ=resultscor.XYZ;   
    
    
     N={'FR2', 'FZ2', 'FCZ', 'CZ1', 'FZA', 'FZ1', 'FL1', 'FL3', 'FL5','CL3', 'CL1', 'CL5', 'PL5', 'PL1', 'PL3','PZC', 'O1', 'PZP', 'OZ', 'O2', 'PR4', 'PR2', 'PR6', 'CR2', 'CZ2', 'CR6', 'CR4', 'FR6', 'FR4'};
       cd(DOI)
    
        STATSFINAL.(textmeasure).go_measure0array=measure0;
        STATSFINAL.(textmeasure).nogo_measure0array=measure0_nogo;
        
        STATSFINAL.(textmeasure).go_measure1array=measure1;
        STATSFINAL.(textmeasure).nogo_measure1array=measure1_nogo;

        STATSFINAL.(textmeasure).go_measure2array=measure2;
        STATSFINAL.(textmeasure).nogo_measure2array=measure2_nogo;

        STATSFINAL.(textmeasure).go_measure3array=measure3_areas;
        STATSFINAL.(textmeasure).nogo_measure3array=measure3_areas_nogo;
            
        STATSFINAL.(textmeasure).go_measure0grandav=grandaverage;
        STATSFINAL.(textmeasure).nogo_measure0granav=grandaverage_nogo;
           
    %excelme2(STATSFINAL, textmeasure, textmeasure, thismoment)
end

% cd(DOI)
% save STATSFINAL STATSFINAL
%% END OF PART 1, CREATION OF STATSFINAL IS DONE.
%% 
% plot2dhead_frontal(grandaverage, XYZ, N); title(['Grandaverage Go-' textmeasure]);
%     figure_temp=['Lines-' textmeasure '- ' 'Grand_average']; 
%     saveas(gcf, figure_temp, 'fig')
%     
% plot2dhead_frontal(grandaverage_nogo, XYZ, N); title(['Grandaverage NoGo-' textmeasure]);
%     figure_temp=['Lines-' textmeasure '- ' 'Grand_averageNoGo']; 
%     saveas(gcf, figure_temp, 'fig')

%     figure; plot2deeg3_frontal(grandaverage, XYZ, N); title(['Grandaverage Go-' textmeasure]);
%     figure_temp=['Grandaverage Go-' textmeasure];
%     saveas(gcf, figure_temp, 'fig');
%     
%     figure; plot2deeg3_frontal(grandaverage_nogo, XYZ,N); title(['Grandaverage NoGo-' textmeasure]);
%     figure_temp=['Grandaverage NoGo-' textmeasure];
%     saveas(gcf, figure_temp, 'fig');
% 
%     for jj=3:18, [p(jj), h(jj)]=ttest2(measure1(jj,:), measure1_nogo(jj,:)); end  %% file-wise statistics
%     figure; plot(h, '*');


% for kk=1:nchan, [hc(kk), pc(kk)]=ttest2(measure1(3:end,kk), measure1_nogo(3:end,kk), 0.005, 'both', 'unequal'); end %% channel-wise statistics
%     figure; plot(pc, '*'); 
%% PART TWO, START OF SOMETHING NEW. 
% GIA KAPOIA TEXTMEASURE - NA GINEI LOOP!!! 23-12-2012
triggerlist={'go', 'nogo'};
% GO
trigger=triggerlist{1}
s=N;
load resultsALL
XYZ=resultsALL.XYZ;
mcor=resultsALL.(textmeasure).go;
numRight=[1,21,22,23,24,26,27,28,29]; % 22-2-2012\\\ 
numLeft=[7:15]; 
% s=resultsALL.s;

nchan=length(XYZ);
corgof=zeros(nchan);
corgofL=zeros(nchan);
corgoall=zeros(nchan);
corgof(numRight, numRight)=mcor(numRight,numRight);
corgofL(numLeft, numLeft)=mcor(numLeft,numLeft); 
corgoall=(corgof+corgofL)./2;
thr_go=0.5*(max(squeeze(max(corgoall))));% to miso tou megistou
%titi=1-0.05^(1/(Ntrial_go-1)); % as Haliday 1995


%% for NOGO 
trigger=triggerlist{2};
s=N;
XYZ=resultsALL.XYZ;
%textmeasure='DTFtheta';
mcor_nogo=resultsALL.(textmeasure).nogo;

nchan=length(XYZ);
corgof_nogo=zeros(nchan);
corgofL_nogo=zeros(nchan);
corgoall_nogo=zeros(nchan);
corgof_nogo(numRight, numRight)=mcor_nogo(numRight,numRight);
corgofL_nogo(numLeft, numLeft)=mcor_nogo(numLeft,numLeft); 
corgoall_nogo=(corgof_nogo+corgofL_nogo)./2;
thr_nogo=0.8*(max(squeeze(max(corgoall_nogo))));% to 80% tou megistou
%titi=1-0.05^(1/(Ntrial_nogo-1)); % as Haliday 1995

thr=min(thr_nogo, thr_go)
grand_lines(textmeasure,corgof,corgofL, XYZ, s, DOI, thr_go, numLeft, numRight,triggerlist{1});
grand_lines(textmeasure,corgof_nogo,corgofL_nogo, XYZ, s, DOI, thr, numLeft, numRight,triggerlist{2});
    
% edw 23-2-2012
%% Here we extract the coupling between brain areas for the grandaverage
cd(DOI)
for k=3:length(textmeasuresall)
    textmeasure=textmeasuresall{k};
    Go1=STATSFINAL.(textmeasure).go_measure0grandav;
    NoGo1=STATSFINAL.(textmeasure).nogo_measure0granav;
       
    Go1_left=Go1(numLeft, numLeft);
    Go1_right=Go1(numRight, numRight);
    NoGo1_left=NoGo1(numLeft, numLeft);
    NoGo1_right=NoGo1(numRight, numRight);
    
    % Left hemisphere 
    %Nleft=N(numLeft)
    Brain_areas={'FL', 'PL', 'CL'}%, 'CZ', 'OZ', 'FZ'};% {'FL', 'PL', 'CL', 'CZ', 'OZ', 'FZ'}
    FL=[1, 2, 3]; PL=[7,8,9]; CL=[4,5,6]; %noumera -indexes- sto yposynolo 
    Go1_left_areas=makeareas_conn(Go1_left, Brain_areas, FL,PL, CL, textmeasure, 'Go', 'Left');% right_area=makeareas_conn(Go1_right, Brain_areas, FR,PR, CR, CZ, OZ, FZ)
    NoGo1_left_areas=makeareas_conn(NoGo1_left, Brain_areas, FL,PL, CL, textmeasure, 'NoGo', 'Left');
    
    % Right Hemisphere
    %Nright=N(numRight)
    FR=[1,8,9];
    PR=[2,3,4];
    CR=[5,6,7];
%     CZ=8;
%     OZ=3;
%     FZ=2;
    Brain_areas={'FR', 'PR', 'CR'}%, 'CZ', 'OZ', 'FZ'};% 
    Go1_right_areas=makeareas_conn(Go1_right, Brain_areas, FR,PR, CR, textmeasure, 'Go', 'Right');% right_area=makeareas_conn(Go1_right, Brain_areas, FR,PR, CR, CZ, OZ, FZ)
    NoGo1_right_areas=makeareas_conn(NoGo1_right, Brain_areas, FR,PR, CR, textmeasure, 'NoGo', 'Right');
    
    % here 23-2-2012 to continueee and finalize tomorrow!!
    
    STATSFINAL.(textmeasure).Go_left_grandav=Go1_left; % maybe we do not need that variable
     STATSFINAL.(textmeasure).Go_left_areas_grandav=Go1_left_areas;
    STATSFINAL.(textmeasure).Go_right_grandav=Go1_right; % maybe we do not need that variable
     STATSFINAL.(textmeasure).Go_right_areas_grandav=Go1_right_areas;
    STATSFINAL.(textmeasure).NoGo_left_grandav=NoGo1_left; % maybe we do not need that variable
     STATSFINAL.(textmeasure).NoGo_left_areas_grandav=NoGo1_left_areas;
    STATSFINAL.(textmeasure).NoGo_right_grandav=NoGo1_right; % maybe we do not need that variable
     STATSFINAL.(textmeasure).NoGo_right_areas_grandav=NoGo1_right_areas;
end

%% asymmetry index Left - Right
for k=3:length(textmeasuresall)
    textmeasure=textmeasuresall{k}
    Left=STATSFINAL.(textmeasure).Go_left_areas_grandav;
    suma1=sum(sum(Left));
    Right=STATSFINAL.(textmeasure).Go_right_areas_grandav;
    suma2=sum(sum(Right));
    LRA=(suma1-suma2)/(suma1+suma2)
    STATSFINAL.(textmeasure).goLRA=LRA;
end

%% assymetry index Anterior posterior
    Brain_areas={'Ant', 'Pos'};% {'FL', 'PL', 'CL', 'CZ', 'OZ', 'FZ'}
    numAnt=[1,2,3,5:9, 28,29];
    numPos=[13:23];
antpostmatrix=makeareas_connAP(Go1, Brain_areas, numAnt,numPos, textmeasure, trigger)
%% SOS  Here we extract the coupling between brain areas for the non-grandaverage- each file
cd(DOI)
for k=1:length(textmeasuresall)
    textmeasure=textmeasuresall{k};
    Go1=STATSFINAL.(textmeasure).go_measure0array; %13x29x29
    NoGo1=STATSFINAL.(textmeasure).nogo_measure0array;
%     numRight=[1,2,20:29]; % 12 el.  Cental electrodes excluded
%     numLeft=[4, 6:15,17]; % 12 el.
    
    Go1_left=Go1(1:end,numLeft, numLeft);
    Go1_right=Go1(1:end,numRight, numRight);
    NoGo1_left=NoGo1(1:end,numLeft, numLeft);
    NoGo1_right=NoGo1(1:end,numRight, numRight);
    
    % Left hemisphere 
    %Nleft=N(numLeft)
    Brain_areas={'FL', 'PL', 'CL', 'CZ', 'OZ', 'FZ'};% {'FL', 'PL', 'CL', 'CZ', 'OZ', 'FZ'}
    FL=[3, 4, 5]; PL=[9,10,11]; CL=[6,7,8]; CZ=1; FZ=2; OZ=10;
    Go1_left_areas=zeros(13,6,6);NoGo1_left_areas=zeros(13,6,6);
    for kk=1:(ND-2) %number of files 
        temp1=squeeze(Go1_left(kk,:,:));
        Go1_left_areas(kk,:,:)=makeareas_conn(temp1, Brain_areas, FL,PL, CL, CZ, OZ, FZ, textmeasure, 'Go', 'Left');% right_area=makeareas_conn(Go1_right, Brain_areas, FR,PR, CR, CZ, OZ, FZ)
        temp2=squeeze(NoGo1_left(kk,:,:));
        NoGo1_left_areas(kk,:,:)=makeareas_conn(temp2, Brain_areas, FL,PL, CL, CZ, OZ, FZ, textmeasure, 'NoGo', 'Left');
        clear temp1 temp2
    end
    
   for lp=1:6 %number of areas
        for lo=1:6
            temp1=Go1_left_areas(:,lp,lo);
            temp2=NoGo1_left_areas(:,lp,lo);
            [h(lp,lo),p(lp,lo)]=ttest2(temp1, temp2, 0.0005, 'right');clear temp1 temp2
        end
    end
   STATSFINAL.(textmeasure).Go_left_areas_p=p;
   figure; imagesc(p<0.05);axis xy; axis tight;set(gca, 'YTickLabel', Brain_areas); set(gca, 'XTickLabel', Brain_areas);
      set(gca,'Ytick', 1:length(Brain_areas)); set(gca, 'XTick', 1:length(Brain_areas));title([textmeasure '-p value-left']); 
   clear p h lo lp
    
    
    % Right Hemisphere
    %Nright=N(numRight)
    FR=[1,11,12];
    PR=[4,5,6];
    CR=[7,9,10];
    CZ=8;
    OZ=3;
    FZ=2;
    Brain_areas={'FR', 'PR', 'CR', 'CZ', 'OZ', 'FZ'};% why not Pz? 28-12-2011
    Go1_right_areas=zeros(13,6,6);NoGo1_right_areas=zeros(13,6,6);
      for kk=1:13
        temp1=squeeze(Go1_right(kk,:,:));
        Go1_right_areas(kk,:,:)=makeareas_conn(temp1, Brain_areas, FR,PR, CR, CZ, OZ, FZ, textmeasure, 'Go', 'Left');% right_area=makeareas_conn(Go1_right, Brain_areas, FR,PR, CR, CZ, OZ, FZ)
        temp2=squeeze(NoGo1_right(kk,:,:));%% to rerun because i had FL, PL, etc
        NoGo1_right_areas(kk,:,:)=makeareas_conn(temp2, Brain_areas, FR,PR, CR, CZ, OZ, FZ, textmeasure, 'NoGo', 'Left');
        clear temp1 temp2
    end
    

       for lp=1:6
        for lo=1:6
            temp1=Go1_right_areas(:,lp,lo); 
            temp2=NoGo1_right_areas(:,lp,lo);
            [h(lp,lo),p(lp,lo)]=ttest2(temp1, temp2, 0.0005, 'right')
        end
       end
    close all
     STATSFINAL.(textmeasure).Go_right_areas_p=p;
     Brain_areas={'FR', 'PR', 'CR', 'CZ', 'OZ', 'FZ'};
    figure; imagesc(p<0.05);axis xy; axis tight;set(gca, 'YTickLabel', Brain_areas); set(gca, 'XTickLabel', Brain_areas);
      set(gca,'Ytick', 1:length(Brain_areas)); set(gca, 'XTick', 1:length(Brain_areas));title([textmeasure '-p value-right']); 
    
    STATSFINAL.(textmeasure).Go_left_all=Go1_left;
     STATSFINAL.(textmeasure).Go_left_areas=Go1_left_areas;
    STATSFINAL.(textmeasure).Go_right_all=Go1_right;
     STATSFINAL.(textmeasure).Go_right_areas=Go1_right_areas;
    STATSFINAL.(textmeasure).NoGo_left_all=NoGo1_left;
     STATSFINAL.(textmeasure).NoGo_left_areas=NoGo1_left_areas;
    STATSFINAL.(textmeasure).NoGo_right_all=NoGo1_right;
     STATSFINAL.(textmeasure).NoGo_right_areas=NoGo1_right_areas;
     
%      stempp=['STATSFINAL-' textmeasure]; % do not delete this!!!!
%     xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
%     titles={name};
%     xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
%     xlswrite(stempp,N, 'Sheet1', 'B2');
%     xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
%     xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
%     xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
%     xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
     
     
end



 cd(DOI)
    
 
 for kk=1:10, [hc3(kk), pc3(kk)]=ttest2(STATSFINAL.(textmeasure).go_measure3array(:,kk), STATSFINAL.(textmeasure).nogo_measure3array(:, kk),0.0005, 'both', 'unequal'); end %% channel-wise statistics    figure; plot(pc3<0.05, '*'); 
    resultstats(k,:)=pc3;



save STATSFINAL STATSFINAL % panta sto telos
%
% cd(DOI)
%    [hc2, pc2]=ttest2(measure1, measure1_nogo, 0.005, 'both', 'unequal');
%    hist(measure1)
%    hist(measure1_nogo)
%    %%
%   % Send the results to excel file 

% stempp=[textmeasure 'SimpleStatistics4-CSD']; % do not delete this!!!!
% xlswrite(stempp, name, 'Sheet1', 'A1:A1');
% titles={'GO', 'NOGO'};
% xlswrite(stempp, (titles(1)), 'Sheet1', 'A2:A2');
% xlswrite(stempp,s, 'Sheet1', 'B2');
% xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
% xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');