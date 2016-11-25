for kkm=1:length(textmeasuresall)
     %% brain areas - INTRAHEMISHPHERIC-AREAS. Here we extract the coupling between brain AREAS for the grandaverage
    textmeasure=textmeasuresall{kkm}
     cd(DOI)
    name='Grandaverage';
    Go1=STATSFINAL(ND-1).(textmeasure).go_measure0grandav;
    NoGo1=STATSFINAL(ND-1).(textmeasure).nogo_measure0granav;
    thr=STATSFINAL(ND-1).(textmeasure).thr;
    
    Go1_left=Go1(numLeft, numLeft);
    Go1_right=Go1(numRight, numRight);
    NoGo1_left=NoGo1(numLeft, numLeft);
    NoGo1_right=NoGo1(numRight, numRight);
    
        
    % Left hemisphere 
    side='Left'; trigger='Go';  Brain_areas={'FL','CL', 'PL'};%
    FL=[1, 2, 3]; PL=[7,8,9]; CL=[4,5,6]; %noumera -indexes- sto yposynolo 
    Go1_left_areas=makeareas_conn(Go1_left, Brain_areas, FL, CL, PL, textmeasure, 'Go', name, DOI, 'Left');
    
    trigger='NoGo';
    NoGo1_left_areas=makeareas_conn(NoGo1_left, Brain_areas, FL,CL, PL, textmeasure, 'NoGo',name, DOI, 'Left');
    
    % Right Hemisphere
    side='Right'; trigger='Go'; FR=[1,8,9]; PR=[2,3,4]; CR=[5,6,7];
    Brain_areas={'FR','CR','PR'}; %
    Go1_right_areas=makeareas_conn(Go1_right, Brain_areas, FR, CR, PR, textmeasure, 'Go', name, DOI, 'Right');
    
    trigger='NoGo';
    NoGo1_right_areas=makeareas_conn(NoGo1_right, Brain_areas, FR, CR, PR, textmeasure, 'NoGo', name, DOI,'Right');
          
    % new threshold for the areas. 
    go_thr_left=max(max(Go1_left_areas))
    go_thr_right=max(max(Go1_right_areas))
    thr_go=(1/2)*min(go_thr_right, go_thr_left)
    
    nogo_thr_left=max(max(NoGo1_left_areas));
    nogo_thr_right=max(max(NoGo1_right_areas));
    thr_nogo=(1/2)*min(nogo_thr_left, nogo_thr_right);
    
    thr=min(thr_go,thr_nogo)
    %% plot again with same threshold
    cd(DOI)
    mkdir('Pictures2')
    cd('Pictures2')
    % Left hemisphere GO
    trigger='Go'; side='Left'; Brain_areas={'FL','CL', 'PL'};%
    FL=[1, 2, 3]; PL=[7,8,9]; CL=[4,5,6]; %noumera -indexes- sto yposynolo 
    figtemp=['Br-Areas-Lines-' name '-' textmeasure '-' trigger '-' side];
    figure; plot3D_frontal(Go1_left, XYZ, numLeft, FL, CL, PL, thr); hold on;
%     title(figtemp);
%     saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
    
        % Right Hemisphere GO
    trigger='Go'; side='Right'; FR=[1,8,9]; PR=[2,3,4]; CR=[5,6,7];
    Brain_areas={'FR','CR','PR'}; %
    figtemp=['Br-Areas-Lines-' name '-' textmeasure '-' trigger];
    plot3D_frontal(Go1_right, XYZ, numRight,FR, CR, PR, thr); title(figtemp);
    saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
    
    
    
    % Left hemisphere nOGO
    trigger='NoGo'; side='Left'; Brain_areas={'FL','CL', 'PL'};%
    figtemp=['Br-Areas-Lines-' name '-' textmeasure '-' trigger '-' side];
    figure; plot3D_frontal(NoGo1_left, XYZ, numLeft,FL, CL, PL, thr); hold on;
%     title(figtemp);
%     saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
   

    % Right Hemisphere nOGO
    trigger='NoGo';side='Right'; FR=[1,8,9]; PR=[2,3,4]; CR=[5,6,7];
    figtemp=['Br-Areas-Lines-' name '-' textmeasure '-' trigger];
    plot3D_frontal(NoGo1_right, XYZ, numRight,FR, CR, PR, thr); title(figtemp)
    saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
    
    cd(DOI)
    STATSFINAL(ND-1).(textmeasure).Go_left_grandav=Go1_left; % maybe we do not need that variable
     STATSFINAL(ND-1).(textmeasure).Go_left_areas_grandav=Go1_left_areas;
    STATSFINAL(ND-1).(textmeasure).Go_right_grandav=Go1_right; % maybe we do not need that variable
     STATSFINAL(ND-1).(textmeasure).Go_right_areas_grandav=Go1_right_areas;
    STATSFINAL(ND-1).(textmeasure).NoGo_left_grandav=NoGo1_left; % maybe we do not need that variable
     STATSFINAL(ND-1).(textmeasure).NoGo_left_areas_grandav=NoGo1_left_areas;
    STATSFINAL(ND-1).(textmeasure).NoGo_right_grandav=NoGo1_right; % maybe we do not need that variable
     STATSFINAL(ND-1).(textmeasure).NoGo_right_areas_grandav=NoGo1_right_areas;
     STATSFINAL(ND-1).(textmeasure).thr_br_areas=thr;
     close all
     clear Go1_left Go1_left_areas Go1_right Go1_right_areas NoGo1_left NoGo1_left_areas NoGo1_right NoGo1_right_areas
end