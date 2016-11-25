function excelme_frontal(STATSFINAL, Sheet1, textmeasure, thismoment,s)% 
%Send the results to excel file 
 
stempp=[thismoment 'stats-grandaverage' textmeasure]; % do not delete this!!!!
    xlswrite(stempp, {thismoment}, 'Sheet1', 'A1:A1');
    titles={'measure1', 'measure2', 'measure 3'};
    xlswrite(stempp, (titles(1)), 'Sheet1', 'A2:A2');
    xlswrite(stempp, {'Go'}, 'Sheet1', 'A3:A3');
    xlswrite(stempp, {s(:)'}, 'Sheet1', 'B3');%% edw den grafei ta onomata twn ilektrodiwn
    xlswrite(stempp,STATSFINAL(14).(textmeasure).go_measure1array,'Sheet1', 'B4');
    N=size(STATSFINAL(14).(textmeasure).go_measure1array,1);
    numtemp1=4+N+2+1; % 4 comes from B4
    numtempf=['B' num2str(numtemp1)];  % B12
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
    xlswrite(stempp, STATSFINAL(14).(textmeasure).nogo_measure1array, 'Sheet1', numtempf);
    xlswrite(stempp,{'NoGo'},'Sheet1', numtempA)
    % measure 3
    numtemp1=numtemp1+N+2+1;  % edw
    numtempf=['B' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp,STATSFINAL(14).(textmeasure).go_measure3array,'Sheet1', numtempf);
    xlswrite(stempp,{'Go Brain-all-areas'},'Sheet1', numtempA)
    
    numtemp1=numtemp1+N+2+1;
    numtempf=['B' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp,STATSFINAL(14).(textmeasure).nogo_measure3array, 'Sheet1', numtempf);
    xlswrite(stempp,{'NoGo Brain-all-areas'},'Sheet1', numtempA)
    
    % Left areas
    % GO
    numtemp1=numtemp1+N+2+1;
    Brain_areas={'FL','CL','PL'};
    numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp,{'Go Left grandaverage 3 areas'},'Sheet1', numtempA)
    xlswrite(stempp, (Brain_areas), 'Sheet1', numtempC);
    
    %arrays
    numtemp1=numtemp1+1;
    numtempB=['B' num2str(numtemp1)]; 
    numtempf=['C' num2str(numtemp1)];
    xlswrite(stempp,(Brain_areas)','Sheet1', numtempB)
    xlswrite(stempp,STATSFINAL(14).(textmeasure).Go_left_areas_grandav, 'Sheet1', numtempf);
   
    %NoGO
     N=size(STATSFINAL(14).(textmeasure).Go_left_areas_grandav,1); % 3
     numtemp1=numtemp1+N+2+1;
    %titles
    numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp,{'NoGo Left grandaverage 3 areas'},'Sheet1', numtempA)
    xlswrite(stempp, (Brain_areas), 'Sheet1', numtempC);
    
    %arrays
    numtemp1=numtemp1+1;
    numtempB=['B' num2str(numtemp1)]; 
    numtempf=['C' num2str(numtemp1)];
    xlswrite(stempp,(Brain_areas)','Sheet1', numtempB)
    xlswrite(stempp,STATSFINAL(14).(textmeasure).NoGo_left_areas_grandav, 'Sheet1', numtempf);
    % here 25-02-2012
    
    % RIGHT HEMISPHERE
    % Left areas
    % GO
    numtemp1=numtemp1+N+2+1;
    Brain_areas={'FR','CR','PR'};
    numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp,{'Go Right grandaverage 3 areas'},'Sheet1', numtempA)
    xlswrite(stempp, (Brain_areas), 'Sheet1', numtempC);
    
    %arrays
    numtemp1=numtemp1+1;
    numtempB=['B' num2str(numtemp1)]; 
    numtempf=['C' num2str(numtemp1)];
    xlswrite(stempp,(Brain_areas)','Sheet1', numtempB)
    xlswrite(stempp,STATSFINAL(14).(textmeasure).Go_right_areas_grandav, 'Sheet1', numtempf);
   
    %NoGO
     N=size(STATSFINAL(14).(textmeasure).NoGo_right_areas_grandav,1); % 3
     numtemp1=numtemp1+N+2+1;
    %titles
    numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp,{'NoGo Right grandaverage 3 areas'},'Sheet1', numtempA)
    xlswrite(stempp, (Brain_areas), 'Sheet1', numtempC);
    
    %arrays
    numtemp1=numtemp1+1;
    numtempB=['B' num2str(numtemp1)]; 
    numtempf=['C' num2str(numtemp1)];
    xlswrite(stempp,(Brain_areas)','Sheet1', numtempB)
    xlswrite(stempp,STATSFINAL(14).(textmeasure).NoGo_right_areas_grandav, 'Sheet1', numtempf);
    
    ND=size(STATSFINAL(14).(textmeasure).go_measure0array,1);
    %% loop for every file
    for kkm=1:(ND-2)
        name=STATSFINAL(kkm).(textmeasure).name(1:5);
        stempp=[name '-' thismoment 'stats' textmeasure]; % do not delete this!!!!
        xlswrite(stempp, {thismoment}, 'Sheet1', 'A1:A1');
        titles={'measure1', 'measure2', 'measure 3'};
        xlswrite(stempp, (titles(1)), 'Sheet1', 'A2:A2');
        
        % Left areas
        % GO
        numtemp1=numtemp1+N+2+1;
        Brain_areas={'FL','CL','PL'};
        numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
        numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
        xlswrite(stempp,{'Go Left grandaverage 3 areas'},'Sheet1', numtempA)
        xlswrite(stempp, (Brain_areas), 'Sheet1', numtempC);

        %arrays
        numtemp1=numtemp1+1;
        numtempB=['B' num2str(numtemp1)]; 
        numtempf=['C' num2str(numtemp1)];
        xlswrite(stempp,(Brain_areas)','Sheet1', numtempB)
        xlswrite(stempp,STATSFINAL(kkm).(textmeasure).Go_left_areas, 'Sheet1', numtempf);
   
        %NoGO
         N=size(STATSFINAL(kkm).(textmeasure).Go_left_areas,1); % 3
         numtemp1=numtemp1+N+2+1;
        %titles
        numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
        numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
        xlswrite(stempp,{'NoGo Left grandaverage 3 areas'},'Sheet1', numtempA)
        xlswrite(stempp, (Brain_areas), 'Sheet1', numtempC);
    
        %arrays
        numtemp1=numtemp1+1;
        numtempB=['B' num2str(numtemp1)]; 
        numtempf=['C' num2str(numtemp1)];
        xlswrite(stempp,(Brain_areas)','Sheet1', numtempB)
        xlswrite(stempp,STATSFINAL(kkm).(textmeasure).NoGo_left_areas, 'Sheet1', numtempf);
        % here 25-02-2012
    
        % RIGHT HEMISPHERE
        % Left areas
        % GO
        numtemp1=numtemp1+N+2+1;
        Brain_areas={'FR','CR','PR'};
        numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
        numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
        xlswrite(stempp,{'Go Right grandaverage 3 areas'},'Sheet1', numtempA)
        xlswrite(stempp, (Brain_areas), 'Sheet1', numtempC);

        %arrays
        numtemp1=numtemp1+1;
        numtempB=['B' num2str(numtemp1)]; 
        numtempf=['C' num2str(numtemp1)];
        xlswrite(stempp,(Brain_areas)','Sheet1', numtempB)
        xlswrite(stempp,STATSFINAL(kkm).(textmeasure).Go_right_areas, 'Sheet1', numtempf);

        %NoGO
         N=size(STATSFINAL(kkm).(textmeasure).NoGo_right_areas,1); % 3
         numtemp1=numtemp1+N+2+1;
        %titles
        numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
        numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
        xlswrite(stempp,{'NoGo Right grandaverage 3 areas'},'Sheet1', numtempA)
        xlswrite(stempp, (Brain_areas), 'Sheet1', numtempC);

        %arrays
        numtemp1=numtemp1+1;
        numtempB=['B' num2str(numtemp1)]; 
        numtempf=['C' num2str(numtemp1)];
        xlswrite(stempp,(Brain_areas)','Sheet1', numtempB)
        xlswrite(stempp,STATSFINAL(kkm).(textmeasure).NoGo_right_areas, 'Sheet1', numtempf);
        end
    
    % here 25-02-2012
% %% measure2
% % titles - headings
%     numtemp1=numtemp1+N+2;
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
%     numtempf=['B' num2str(numtemp1) ':B' num2str(numtemp1)]; 
%     numtempC=['C' num2str(numtemp1) ':C' num2str(numtemp1)];  
%     numtempD=['D' num2str(numtemp1) ':D' num2str(numtemp1)];  
%     numtempE=['E' num2str(numtemp1) ':E' num2str(numtemp1)];  
%     xlswrite(stempp, (titles(2)), 'Sheet1', numtempA);
%     xlswrite(stempp, {'Str Sleep'}, 'Sheet1', numtempf);
%     xlswrite(stempp, {'Str Awake'}, 'Sheet1', numtempC);
%     xlswrite(stempp, {'App Awake'}, 'Sheet1', numtempD);
%     xlswrite(stempp, {'App Awake'}, 'Sheet1', numtempE);
% 
% % array
%     numtemp1=numtemp1+1;
%     numtempf=['B' num2str(numtemp1)];
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
%     xlswrite(stempp,STATSFINAL(11).(textmeasure).array2,'Sheet1', numtempf);
% 
% % p value
%     N=length(STATSFINAL(11).(textmeasure).array2);
%     numtemp1=numtemp1+N+1;
%     numtempf=['B' num2str(numtemp1)];  
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
%     xlswrite(stempp,{'p value str'},'Sheet1', numtempA);
%     xlswrite(stempp,STATSFINAL(11).(textmeasure).measure2p_str, 'Sheet1', numtempf);
%     numtemp1=numtemp1+1;
%     numtempf=['B' num2str(numtemp1)];  
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
%     xlswrite(stempp,{'p value app'},'Sheet1', numtempA)
%     xlswrite(stempp,STATSFINAL(11).(textmeasure).measure2p_app, 'Sheet1', numtempf);
% %
% 
% %measure2median1
%     numtemp1=numtemp1+1;
%     numtempf=['B' num2str(numtemp1)];  
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
%     xlswrite(stempp,{'median Str-Sleep'},'Sheet1', numtempA)
%     xlswrite(stempp,STATSFINAL(11).(textmeasure).measure2median1, 'Sheet1', numtempf);
%     numtemp1=numtemp1+1;
%     numtempf=['B' num2str(numtemp1)];  
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
%     xlswrite(stempp,{'Percent Str-Sleep'},'Sheet1', numtempA)
%     xlswrite(stempp,STATSFINAL(11).(textmeasure).measure2percentile1, 'Sheet1', numtempf);
% %measure2median2
%     numtemp1=numtemp1+1;
%     numtempf=['B' num2str(numtemp1)];  
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
%     xlswrite(stempp,{'median Str-Awake'},'Sheet1', numtempA)
%     xlswrite(stempp,STATSFINAL(11).(textmeasure).measure2median2, 'Sheet1', numtempf);
%     numtemp1=numtemp1+1;
%     numtempf=['B' num2str(numtemp1)];  
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
%     xlswrite(stempp,{'Percent Str-Awake'},'Sheet1', numtempA)
%     xlswrite(stempp,STATSFINAL(11).(textmeasure).measure2percentile2, 'Sheet1', numtempf);
% %measure2median3
%     numtemp1=numtemp1+1;
%     numtempf=['B' num2str(numtemp1)];  
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
%     xlswrite(stempp,{'median App-Sleep'},'Sheet1', numtempA)
%     xlswrite(stempp,STATSFINAL(11).(textmeasure).measure2median3, 'Sheet1', numtempf);
%     numtemp1=numtemp1+1;
%     numtempf=['B' num2str(numtemp1)];  
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
%     xlswrite(stempp,{'Percent App-Sleep'},'Sheet1', numtempA)
%     xlswrite(stempp,STATSFINAL(11).(textmeasure).measure2percentile3, 'Sheet1', numtempf);
% %measure2median4
%     numtemp1=numtemp1+1;
%     numtempf=['B' num2str(numtemp1)];  
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
%     xlswrite(stempp,{'median App-Awake'},'Sheet1', numtempA)
%     xlswrite(stempp,STATSFINAL(11).(textmeasure).measure2median4, 'Sheet1', numtempf);
%     numtemp1=numtemp1+1;
%     numtempf=['B' num2str(numtemp1)];  
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
%     xlswrite(stempp,{'Percent App-Awake'},'Sheet1', numtempA)
%     xlswrite(stempp,STATSFINAL(11).(textmeasure).measure2percentile4, 'Sheet1', numtempf);
% 
% %% measure 3
% % titles - headings
%     numtemp1=numtemp1+2;
%     numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
%     numtempf=['B' num2str(numtemp1) ':B' num2str(numtemp1)]; 
%     numtempC=['C' num2str(numtemp1) ':C' num2str(numtemp1)];  
%     numtempD=['D' num2str(numtemp1) ':D' num2str(numtemp1)];  
%     numtempE=['E' num2str(numtemp1) ':E' num2str(numtemp1)];  
%     xlswrite(stempp, (titles(3)), 'Sheet1', numtempA);
%     xlswrite(stempp, {'Day'}, 'Sheet1', numtempf);
%     xlswrite(stempp, {'Night'}, 'Sheet1', numtempC);
% xlswrite(stempp, {'Night'}, 'Sheet1', numtempD);
% xlswrite(stempp, {'Day'}, 'Sheet1', numtempE);
% 
% % array
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,STATSFINAL(11).(textmeasure).array3,'Sheet1', numtempf);
% %
% % p value
% N=length(STATSFINAL(11).(textmeasure).array3);
% numtemp1=numtemp1+N+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
% xlswrite(stempp,{'p value Day2Night'},'Sheet1', numtempA);
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure3p_day2nigth, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
% xlswrite(stempp,{'p value Night2Day'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure3p_nigth2day, 'Sheet1', numtempf);
% 
% % median values & percentiles
% %measure3 median 1
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median Day'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure3p_1_median, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent Day'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure3p_1_percent, 'Sheet1', numtempf);
% %measure3 median 2
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median Night'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure3p_2_median, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent Night'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure3p_2_percent, 'Sheet1', numtempf);
% %measure3 median 3
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median Night'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure3p_3_median, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent Night'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure3p_3_percent, 'Sheet1', numtempf);
% %measure3 median 4
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median Day'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure3p_4_median, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent Day'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure3p_4_percent, 'Sheet1', numtempf);
% 
% %% measure 4
% % titles - headings
% numtemp1=numtemp1+2;
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% numtempf=['B' num2str(numtemp1) ':B' num2str(numtemp1)]; 
% numtempC=['C' num2str(numtemp1) ':C' num2str(numtemp1)];  
% numtempD=['D' num2str(numtemp1) ':D' num2str(numtemp1)];  
% numtempE=['E' num2str(numtemp1) ':E' num2str(numtemp1)];  
% numtempF=['F' num2str(numtemp1) ':F' num2str(numtemp1)]; 
% numtempG=['G' num2str(numtemp1) ':G' num2str(numtemp1)];  
% numtempH=['H' num2str(numtemp1) ':H' num2str(numtemp1)];  
% numtempI=['I' num2str(numtemp1) ':I' num2str(numtemp1)];  
% xlswrite(stempp, (titles(4)), 'Sheet1', numtempA);
% xlswrite(stempp, {'str sleep spik'}, 'Sheet1', numtempf);
% xlswrite(stempp, {'str sleep max'}, 'Sheet1', numtempC);
% xlswrite(stempp, {'str awake spik'}, 'Sheet1', numtempD);
% xlswrite(stempp, {'str awake max'}, 'Sheet1', numtempE);
% xlswrite(stempp, {'app sleep spik'}, 'Sheet1', numtempF);
% xlswrite(stempp, {'app sleep max'}, 'Sheet1', numtempG);
% xlswrite(stempp, {'app awake spik'}, 'Sheet1', numtempH);
% xlswrite(stempp, {'app awake max'}, 'Sheet1', numtempI);
% 
% % array
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,STATSFINAL(11).(textmeasure).array4,'Sheet1', numtempf);
% %
% % p value
% N=length(STATSFINAL(11).(textmeasure).array4);
% numtemp1=numtemp1+N+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
% xlswrite(stempp,{'p value str_max_spike_day'},'Sheet1', numtempA);
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_str_max_spike_day, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
% xlswrite(stempp,{'p value str_max_spike_night'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_str_max_spike_night, 'Sheet1', numtempf);
% % add 2012
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
% xlswrite(stempp,{'p value app_max_spike_day'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_app_max_spike_day, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
% xlswrite(stempp,{'p value app_max_spike_night'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_app_max_spike_night, 'Sheet1', numtempf);
% 
% 
% % median values & percentiles
% %measure4 median 1
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median Str night spik'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_median1, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent Str night spik'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_percent1, 'Sheet1', numtempf);
% %measure4 median 2
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median str night max'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_median2, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent str night max'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_percent2, 'Sheet1', numtempf);
% %measure4 median 3
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median str day spik'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_median3, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent str day spik'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_percent3, 'Sheet1', numtempf);
% %measure4 median 4
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median str day max'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_median4, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent str day max'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_percent4, 'Sheet1', numtempf);
% %measure4 median 5
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median app night spik'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_median5, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent app night spik'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_percent5, 'Sheet1', numtempf);
% %measure4 median 6
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median app night max'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_median6, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent app night max'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_percent6, 'Sheet1', numtempf);
% %measure4 median 7
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median app day spik'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_median7, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent app day spik'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_percent7, 'Sheet1', numtempf);
% %measure4 median 8
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median app day max'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_median8, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent app day max'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure4p_percent8, 'Sheet1', numtempf);
% %% measure 5
% % titles - headings
% numtemp1=numtemp1+2;
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% numtempf=['B' num2str(numtemp1) ':B' num2str(numtemp1)]; 
% numtempC=['C' num2str(numtemp1) ':C' num2str(numtemp1)];  
% numtempD=['D' num2str(numtemp1) ':D' num2str(numtemp1)];  
% numtempE=['E' num2str(numtemp1) ':E' num2str(numtemp1)];  
% xlswrite(stempp, (titles(5)), 'Sheet1', numtempA);
% xlswrite(stempp, {'Str sleep Occ'}, 'Sheet1', numtempf);
% xlswrite(stempp, {'Str awake Occ'}, 'Sheet1', numtempC);
% xlswrite(stempp, {'app sleep Occ'}, 'Sheet1', numtempD);
% xlswrite(stempp, {'app awake Occ'}, 'Sheet1', numtempE);
% 
% % array
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];
% xlswrite(stempp,STATSFINAL(11).(textmeasure).array5,'Sheet1', numtempf);
% 
% %
% % p value
% N=length(STATSFINAL(11).(textmeasure).array5);
% numtemp1=numtemp1+N+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
% xlswrite(stempp,{'p value Occ str '},'Sheet1', numtempA);
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure5p_str, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
% xlswrite(stempp,{'p value Occ app'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure5p_app, 'Sheet1', numtempf);
% 
% % median values & percentiles
% %measure5 median 1
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median Occ sleep str'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure5p_median1, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent Occ sleep str'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure5p_percent1, 'Sheet1', numtempf);
% %measure5 median 2
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median Occ str awake'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure5p_median2, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'percent Occ str awake'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure5p_percent2, 'Sheet1', numtempf);
% %measure5 median 3
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median Occ app sleep'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure5p_median3, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent Occ app sleep'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure5p_percent3, 'Sheet1', numtempf);
% %measure5 median 4
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
% xlswrite(stempp,{'median Opp app awake'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure5p_median4, 'Sheet1', numtempf);
% numtemp1=numtemp1+1;
% numtempf=['B' num2str(numtemp1)];  
% numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
% xlswrite(stempp,{'Percent Opp app awake'},'Sheet1', numtempA)
% xlswrite(stempp,STATSFINAL(11).(textmeasure).measure5p_percent4, 'Sheet1', numtempf);
% 
% % function encore('Sheet1', 'text1', measurexp)
% % numtemp1=numtemp1+1;
% % numtempf=['B' num2str(numtemp1)];  
% % numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
