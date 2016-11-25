% Make statistical analysis 07-12-2011, rev 23-12-2011
clear all 
close all
% cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\Second component')
% DOI='D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\Second component';

%doi='D:\OFC\SETS_correct_locations\'; %'D:\RIKSHOSPITALET\CNV_RIKS\RAW DATASETS\After_eye_art_removal'
%DOI='D:\OFC\ANALYZED_DATASETS\'; 
DOI='D:\OFC\ANALYZED_DATASETS_cube\All_interval_OFC'
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
    % for GO
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
            case {'DTFgamma1', 'DTFdelta'}, chan_appearance_norm=resultscor.(textmeasure).chan_appearance_norm; %meantemp=resultscor.resultsDTF.(textmeasurex).meangamma;
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
        clear chan_strength_norm chan_appearance_norm 
                 
        
        all_areas={'FL', 'FR', 'CZ', 'CL', 'CR', 'PZ','PL','PR', 'OZ', 'FZ'}
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
         
        STATSFINAL(ND-1).(textmeasure).go_measure0array=measure0;
        STATSFINAL(ND-1).(textmeasure).nogo_measure0array=measure0_nogo;
        
        STATSFINAL(ND-1).(textmeasure).go_measure1array=measure1;
        STATSFINAL(ND-1).(textmeasure).nogo_measure1array=measure1_nogo;

        STATSFINAL(ND-1).(textmeasure).go_measure2array=measure2;
        STATSFINAL(ND-1).(textmeasure).nogo_measure2array=measure2_nogo;

        STATSFINAL(ND-1).(textmeasure).go_measure3array=measure3_areas;
        STATSFINAL(ND-1).(textmeasure).nogo_measure3array=measure3_areas_nogo;
            
        STATSFINAL(ND-1).(textmeasure).go_measure0grandav=grandaverage;
        STATSFINAL(ND-1).(textmeasure).nogo_measure0granav=grandaverage_nogo;
        STATSFINAL(ND-1).(textmeasure).all_areas=all_areas;   
    %excelme2(STATSFINAL, textmeasure, textmeasure, thismoment)

    clear resultscor allcouples_nogo allcouples
    clear measure0 measure0_nogo measure1 measure1_nogo measure2 measure2_nogo measure3_areas measure3_areas_nogo grandaverage grandaverage_nogo
%% PART TWO, RIGHT-LEFT CONNECTIVITY INTRA-HEMISPHERIC
% GIA KAPOIA TEXTMEASURE - NA GINEI LOOP!!! 23-12-2012, we can continue the
% above loop with 
% FOR GRANDAVERAGE, AND SHOWING PICTURES PURPOSES
    cd(DOI)
    triggerlist={'go', 'nogo'};
% GO
    trigger=triggerlist{1}
    s=N;
    load resultsALL
    XYZ=resultsALL.XYZ;
    mcor=resultsALL.(textmeasure).go; % edw einai to grandaverage.
    numRight=[1,21,22,23,24,26,27,28,29]; % 22-2-2012\\\ 
    numLeft=[7:15]; 
    % s=resultsALL.s;

    [connL_go connR_go thr_go]=left_right(mcor, numRight, numLeft, XYZ);

%% for NOGO 
    trigger=triggerlist{2};
    s=N;
    mcor_nogo=resultsALL.(textmeasure).nogo;
    [connL_nogo connR_nogo thr_nogo]=left_right(mcor_nogo, numRight, numLeft, XYZ);

    % COMMON THRESHOLD BETWEEN GO AND NOGO ASKED BY ANNE-KRISTIN
    thr=min(thr_nogo, thr_go);
    name='Grandaverage'
    grand_lines(textmeasure,name, connL_go,connR_go, XYZ, s, DOI, thr, numLeft, numRight,triggerlist{1});
    grand_lines(textmeasure,name,connL_nogo,connR_nogo, XYZ, s, DOI, thr, numLeft, numRight,triggerlist{2});
    STATSFINAL(ND-1).(textmeasure).name=name; 
    STATSFINAL(ND-1).(textmeasure).connL_nogo=connL_nogo; 
    STATSFINAL(ND-1).(textmeasure).connR_nogo=connR_nogo; 
    STATSFINAL(ND-1).(textmeasure).connL_go=connL_go; 
    STATSFINAL(ND-1).(textmeasure).connR_go=connR_go; 
    STATSFINAL(ND-1).(textmeasure).thr_go=thr_go; 
    STATSFINAL(ND-1).(textmeasure).thr_nogo=thr_nogo; 
    STATSFINAL(ND-1).(textmeasure).thr=thr; 
    close all
    clear connL_nogo connR_nogo connL connR thr_go thr_nogo thr mcor_nogo mcor_go name kkm
% INTRAHEMISPHERIC for every file
    for kkm=1:(ND-2)
        name=filenames{kkm+2};
        mcor_go(:,:)=STATSFINAL(ND-1).(textmeasure).go_measure0array(kkm,:,:);
        mcor_go=squeeze(mcor_go);
        [connL_go connR_go thr_go]=left_right(mcor_go, numRight, numLeft, XYZ);

        mcor_nogo=STATSFINAL(ND-1).(textmeasure).nogo_measure0array(kkm,:,:);
        mcor_nogo=squeeze(mcor_go);
        [connL_nogo connR_nogo thr_nogo]=left_right(mcor_nogo, numRight, numLeft, XYZ);

        thr=min(thr_nogo, thr_go);
        grand_lines(textmeasure,name, connL_go,connR_go, XYZ, s, DOI, thr, numLeft, numRight, triggerlist{1});
        grand_lines(textmeasure,name, connL_nogo,connR_nogo, XYZ, s, DOI, thr, numLeft, numRight,triggerlist{2});
        
        STATSFINAL(kkm).(textmeasure).name=name; 
        STATSFINAL(kkm).(textmeasure).connL_nogo=connL_nogo; 
        STATSFINAL(kkm).(textmeasure).connR_nogo=connR_nogo; 
        STATSFINAL(kkm).(textmeasure).connL_go=connL_go; 
        STATSFINAL(kkm).(textmeasure).connR_go=connR_go; 
        STATSFINAL(kkm).(textmeasure).thr_go=thr_go; 
        STATSFINAL(kkm).(textmeasure).thr_nogo=thr_nogo; 
        STATSFINAL(kkm).(textmeasure).thr=thr;     
    end
        clear connL_nogo connR_nogo connL_go connR_go thr_go thr_nogo thr mcor_nogo mcor_go mcor
        close all
    % edw 23-2-2012
        clear figtemp a filenames2 files2 kk qq resultsALL
    %% brain areas - INTRAHEMISHPHERIC-AREAS. Here we extract the coupling between brain AREAS for the grandaverage
    Go1=STATSFINAL(ND-1).(textmeasure).go_measure0grandav;
    NoGo1=STATSFINAL(ND-1).(textmeasure).nogo_measure0granav;
    thr=STATSFINAL(ND-1).(textmeasure).thr
    
    Go1_left=Go1(numLeft, numLeft);
    Go1_right=Go1(numRight, numRight);
    NoGo1_left=NoGo1(numLeft, numLeft);
    NoGo1_right=NoGo1(numRight, numRight);
    
    % Left hemisphere 
    side='Left'; trigger='Go';  Brain_areas={'FL','CL', 'PL'};%
    FL=[1, 2, 3]; PL=[7,8,9]; CL=[4,5,6]; %noumera -indexes- sto yposynolo 
    Go1_left_areas=makeareas_conn(Go1_left, Brain_areas, FL, CL, PL, textmeasure, 'Go', name, DOI, 'Left');
    figtemp=['Areas-Lines-' name(1:5) '-' textmeasure '-' trigger '-' side];
    plot3D_frontal(Go1_left, XYZ, numLeft, FL, CL, PL, thr); title(figtemp);
    saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
    
    trigger='NoGo';
    NoGo1_left_areas=makeareas_conn(NoGo1_left, Brain_areas, FL,CL, PL, textmeasure, 'NoGo',name, DOI, 'Left');
    figtemp=['Areas-Lines-' name(1:5) '-' textmeasure '-' trigger '-' side];
    plot3D_frontal(NoGo1_left, XYZ, numLeft,FL, CL, PL, thr); 
    title(figtemp);
    saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
   
    % Right Hemisphere
    side='Right'; trigger='Go'; FR=[1,8,9]; PR=[2,3,4]; CR=[5,6,7];
    Brain_areas={'FR','CR','PR'}; %
    Go1_right_areas=makeareas_conn(Go1_right, Brain_areas, FR, CR, PR, textmeasure, 'Go', name, DOI, 'Right');
    figtemp=['Areas-Lines-' name(1:5) '-' textmeasure '-' trigger '-' side];
    plot3D_frontal(Go1_right, XYZ, numRight,FR, CR, PR, thr); title(figtemp);
    saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
    
    trigger='NoGo';
    NoGo1_right_areas=makeareas_conn(NoGo1_right, Brain_areas, FR, CR, PR, textmeasure, 'NoGo', name, DOI,'Right');
    figtemp=['Areas-Lines-' name(1:5) '-' textmeasure '-' trigger '-' side];
    plot3D_frontal(NoGo1_right, XYZ, numRight,FR, CR, PR, thr); title(figtemp)
    saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
    
    STATSFINAL(ND-1).(textmeasure).Go_left_grandav=Go1_left; % maybe we do not need that variable
     STATSFINAL(ND-1).(textmeasure).Go_left_areas_grandav=Go1_left_areas;
    STATSFINAL(ND-1).(textmeasure).Go_right_grandav=Go1_right; % maybe we do not need that variable
     STATSFINAL(ND-1).(textmeasure).Go_right_areas_grandav=Go1_right_areas;
    STATSFINAL(ND-1).(textmeasure).NoGo_left_grandav=NoGo1_left; % maybe we do not need that variable
     STATSFINAL(ND-1).(textmeasure).NoGo_left_areas_grandav=NoGo1_left_areas;
    STATSFINAL(ND-1).(textmeasure).NoGo_right_grandav=NoGo1_right; % maybe we do not need that variable
     STATSFINAL(ND-1).(textmeasure).NoGo_right_areas_grandav=NoGo1_right_areas;
     close all
     clear Go1_left Go1_left_areas Go1_right Go1_right_areas NoGo1_left NoGo1_left_areas NoGo1_right NoGo1_right_areas
      %% Here we extract the coupling between brain AREAS for each file
      for kkm=1:(ND-2)
          name=filenames{kkm+2};
          Go1(:,:)=STATSFINAL(ND-1).(textmeasure).go_measure0array(kkm,:,:);
          NoGo1(:,:)=STATSFINAL(ND-1).(textmeasure).nogo_measure0array(kkm,:,:);
          thr=STATSFINAL(ND-1).(textmeasure).thr;
          
          Go1_left=Go1(numLeft, numLeft);
          Go1_right=Go1(numRight, numRight);
          NoGo1_left=NoGo1(numLeft, numLeft);
          NoGo1_right=NoGo1(numRight, numRight);
         
          % Left hemisphere 
    	  side='Left'; Brain_areas={'FL', 'CL', 'PL'}; trigger='Go';
          Go1_left_areas=makeareas_conn(Go1_left, Brain_areas, FL, CL, PL, textmeasure, 'Go', name, DOI, 'Left');
          %figtemp=['Areas-Lines-' name '-' textmeasure '-' trigger '-' side];
          % plot3D_frontal(Go1_left, XYZ, numLeft,Far, Car, Par, thr); 
          %saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
    
          trigger='NoGo';
          NoGo1_left_areas=makeareas_conn(NoGo1_left, Brain_areas, FL,CL, PL, textmeasure, 'NoGo',name, DOI, 'Left');
          %figtemp=['Areas-Lines-' name '-' textmeasure '-' trigger '-' side];
          %plot3D_frontal(NoGo1_left, XYZ, numLeft,Far, Car, Par, thr); 
          %saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
   
          % Right Hemisphere
            side='Right'; Brain_areas={'FR','CR', 'PR'}; trigger='Go';
            FR=[1,8,9]; PR=[2,3,4]; CR=[5,6,7];
            %
            Go1_right_areas=makeareas_conn(Go1_right, Brain_areas, FR, CR, PR, textmeasure, 'Go', name, DOI, 'Right');
            %figtemp=['Areas-Lines-' name '-' textmeasure '-' trigger '-' side];
            %plot3D_frontal(Go1_right, XYZ, numRight,Far, Car, Par, thr); 
            %saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
    
            trigger='NoGo';
            NoGo1_right_areas=makeareas_conn(NoGo1_right, Brain_areas, FR, CR, PR, textmeasure, 'NoGo', name, DOI,'Right');
            %figtemp=['Areas-Lines-' name '-' textmeasure '-' trigger '-' side];
            %plot3D_frontal(NoGo1_right, XYZ, numLeft,Far, Car, Par, thr); 
            %saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
    
            STATSFINAL(kkm).(textmeasure).Go_left=Go1_left; % maybe we do not need that variable
            STATSFINAL(kkm).(textmeasure).Go_left_areas=Go1_left_areas;
            STATSFINAL(kkm).(textmeasure).Go_right=Go1_right; % maybe we do not need that variable
            STATSFINAL(kkm).(textmeasure).Go_right_areas=Go1_right_areas;
            STATSFINAL(kkm).(textmeasure).NoGo_left=NoGo1_left; % maybe we do not need that variable
            STATSFINAL(kkm).(textmeasure).NoGo_left_areas=NoGo1_left_areas;
            STATSFINAL(kkm).(textmeasure).NoGo_right=NoGo1_right; % maybe we do not need that variable
            STATSFINAL(kkm).(textmeasure).NoGo_right_areas=NoGo1_right_areas;
            % here 24-2-2012 to continueee and finalize tomorrow!!
      end

    clear Go1_left Go1_left_areas Go1_right Go1_right_areas NoGo1_left_areas NoGo1_right NoGo1_right_areas
    close all
    cd(DOI)
    excelme_frontal2(STATSFINAL, textmeasure textmeasure, thismoment,s)% 
% %% asymmetry index Left - Right for the grandaverage
%     for k=3:length(textmeasuresall)
%         textmeasure=textmeasuresall{k}
%         Left=STATSFINAL.(textmeasure).Go_left_areas_grandav;
%         suma1=sum(sum(Left));
%         Right=STATSFINAL.(textmeasure).Go_right_areas_grandav;
%         suma2=sum(sum(Right));
%         LRA=(suma1-suma2)/(suma1+suma2)
%         STATSFINAL.(textmeasure).goLRA=LRA;
%     end

% %% assymetry index Anterior posterior
%     Brain_areas={'Ant', 'Pos'};% {'FL', 'PL', 'CL', 'CZ', 'OZ', 'FZ'}
%     numAnt=[1,2,3,5:9, 28,29];
%     numPos=[13:23];
% antpostmatrix=makeareas_connAP(Go1, Brain_areas, numAnt,numPos, textmeasure, trigger)
end

save STATSFINAL STATSFINAL

% %% SOS  Here we extract the coupling between brain areas for the non-grandaverage- each file
% cd(DOI)
% for k=1:length(textmeasuresall)
%     textmeasure=textmeasuresall{k};
%     Go1=STATSFINAL.(textmeasure).go_measure0array; %13x29x29
%     NoGo1=STATSFINAL.(textmeasure).nogo_measure0array;
% %     numRight=[1,2,20:29]; % 12 el.  Cental electrodes excluded
% %     numLeft=[4, 6:15,17]; % 12 el.
%     
%     Go1_left=Go1(1:end,numLeft, numLeft);
%     Go1_right=Go1(1:end,numRight, numRight);
%     NoGo1_left=NoGo1(1:end,numLeft, numLeft);
%     NoGo1_right=NoGo1(1:end,numRight, numRight);
%     
%     % Left hemisphere 
%     %Nleft=N(numLeft)
%     Brain_areas={'FL', 'PL', 'CL', 'CZ', 'OZ', 'FZ'};% {'FL', 'PL', 'CL', 'CZ', 'OZ', 'FZ'}
%     FL=[3, 4, 5]; PL=[9,10,11]; CL=[6,7,8]; CZ=1; FZ=2; OZ=10;
%     Go1_left_areas=zeros(13,6,6);NoGo1_left_areas=zeros(13,6,6);
%     for kk=1:(ND-2) %number of files 
%         temp1=squeeze(Go1_left(kk,:,:));
%         Go1_left_areas(kk,:,:)=makeareas_conn(temp1, Brain_areas, FL,PL, CL, CZ, OZ, FZ, textmeasure, 'Go', 'Left');% right_area=makeareas_conn(Go1_right, Brain_areas, FR,PR, CR, CZ, OZ, FZ)
%         temp2=squeeze(NoGo1_left(kk,:,:));
%         NoGo1_left_areas(kk,:,:)=makeareas_conn(temp2, Brain_areas, FL,PL, CL, CZ, OZ, FZ, textmeasure, 'NoGo', 'Left');
%         clear temp1 temp2
%     end
%     
%    for lp=1:6 %number of areas
%         for lo=1:6
%             temp1=Go1_left_areas(:,lp,lo);
%             temp2=NoGo1_left_areas(:,lp,lo);
%             [h(lp,lo),p(lp,lo)]=ttest2(temp1, temp2, 0.0005, 'right');clear temp1 temp2
%         end
%     end
%    STATSFINAL.(textmeasure).Go_left_areas_p=p;
%    figure; imagesc(p<0.05);axis xy; axis tight;set(gca, 'YTickLabel', Brain_areas); set(gca, 'XTickLabel', Brain_areas);
%       set(gca,'Ytick', 1:length(Brain_areas)); set(gca, 'XTick', 1:length(Brain_areas));title([textmeasure '-p value-left']); 
%    clear p h lo lp
%     
%     
%     % Right Hemisphere
%     %Nright=N(numRight)
%     FR=[1,11,12];
%     PR=[4,5,6];
%     CR=[7,9,10];
%     CZ=8;
%     OZ=3;
%     FZ=2;
%     Brain_areas={'FR', 'PR', 'CR', 'CZ', 'OZ', 'FZ'};% why not Pz? 28-12-2011
%     Go1_right_areas=zeros(13,6,6);NoGo1_right_areas=zeros(13,6,6);
%       for kk=1:13
%         temp1=squeeze(Go1_right(kk,:,:));
%         Go1_right_areas(kk,:,:)=makeareas_conn(temp1, Brain_areas, FR,PR, CR, CZ, OZ, FZ, textmeasure, 'Go', 'Left');% right_area=makeareas_conn(Go1_right, Brain_areas, FR,PR, CR, CZ, OZ, FZ)
%         temp2=squeeze(NoGo1_right(kk,:,:));%% to rerun because i had FL, PL, etc
%         NoGo1_right_areas(kk,:,:)=makeareas_conn(temp2, Brain_areas, FR,PR, CR, CZ, OZ, FZ, textmeasure, 'NoGo', 'Left');
%         clear temp1 temp2
%     end
%     
% 
%        for lp=1:6
%         for lo=1:6
%             temp1=Go1_right_areas(:,lp,lo); 
%             temp2=NoGo1_right_areas(:,lp,lo);
%             [h(lp,lo),p(lp,lo)]=ttest2(temp1, temp2, 0.0005, 'right')
%         end
%        end
%     close all
%      STATSFINAL.(textmeasure).Go_right_areas_p=p;
%      Brain_areas={'FR', 'PR', 'CR', 'CZ', 'OZ', 'FZ'};
%     figure; imagesc(p<0.05);axis xy; axis tight;set(gca, 'YTickLabel', Brain_areas); set(gca, 'XTickLabel', Brain_areas);
%       set(gca,'Ytick', 1:length(Brain_areas)); set(gca, 'XTick', 1:length(Brain_areas));title([textmeasure '-p value-right']); 
%     
%     STATSFINAL.(textmeasure).Go_left_all=Go1_left;
%      STATSFINAL.(textmeasure).Go_left_areas=Go1_left_areas;
%     STATSFINAL.(textmeasure).Go_right_all=Go1_right;
%      STATSFINAL.(textmeasure).Go_right_areas=Go1_right_areas;
%     STATSFINAL.(textmeasure).NoGo_left_all=NoGo1_left;
%      STATSFINAL.(textmeasure).NoGo_left_areas=NoGo1_left_areas;
%     STATSFINAL.(textmeasure).NoGo_right_all=NoGo1_right;
%      STATSFINAL.(textmeasure).NoGo_right_areas=NoGo1_right_areas;
%  
%      
% end
% 
% 
% 
%  cd(DOI)
%     
%  
%  for kk=1:10, [hc3(kk), pc3(kk)]=ttest2(STATSFINAL.(textmeasure).go_measure3array(:,kk), STATSFINAL.(textmeasure).nogo_measure3array(:, kk),0.0005, 'both', 'unequal'); end %% channel-wise statistics    figure; plot(pc3<0.05, '*'); 
%     resultstats(k,:)=pc3;
% 
% 
% 
% save STATSFINAL STATSFINAL % panta sto telos
% %
% % cd(DOI)
% %    [hc2, pc2]=ttest2(measure1, measure1_nogo, 0.005, 'both', 'unequal');
% %    hist(measure1)
% %    hist(measure1_nogo)
% %    %%
% %   % Send the results to excel file 
% 
% % stempp=[textmeasure 'SimpleStatistics4-CSD']; % do not delete this!!!!
% % xlswrite(stempp, name, 'Sheet1', 'A1:A1');
% % titles={'GO', 'NOGO'};
% % xlswrite(stempp, (titles(1)), 'Sheet1', 'A2:A2');
% % xlswrite(stempp,s, 'Sheet1', 'B2');
% % xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
% % xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');