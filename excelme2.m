function excelme2(STATSFINAL, Sheet1, textmeasure, thismoment)% 
stempp=[thismoment 'stats' textmeasure]; % do not delete this!!!!
    xlswrite(stempp, {thismoment}, 'Sheet1', 'A1:A1');
    titles={'grandaverage','measure1', 'measure2', 'measure3'};
   %LINESA={'A', 'B', 'C', 'D', 'E'''''''''''''''''''''''}
    % measure 0 
    xlswrite(stempp, (titles(1)), 'Sheet1', 'A2:A2');
    xlswrite(stempp, {'GO'}, 'Sheet1', 'A3:A3');
    xlswrite(stempp, STATSFINAL(14).(textmeasure).go_measure0grandav, 'Sheet1', 'A4');
    N=length(STATSFINAL.(textmeasure).go_measure0grandav);
    numtemp1=N+4+1; % TO 4 OFEILETAI STO GEGONOS OTI EXOUME FTASEI MEXRI TO A4
    numtempf=['B' num2str(numtemp1)];  % B12
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp, {'NOGO'}, 'Sheet1', numtempA);
    xlswrite(stempp,STATSFINAL.go_measure0grandav,'Sheet1', numtempf);
    numtemp1=numtemp1+1+N;
    % measure 1
    numtempf=['B' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp,{'Measure1 channel strength GO'},'Sheet1', numtempA);
    xlswrite(stempp,STATSFINAL.(textmeasure).go_measure1array,'Sheet1', numtempf);
    N=size(STATSFINAL.(textmeasure).go_measure1array,1)
    numtemp1=numtemp1+N+1;
    
    
    numtempf=['B' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp,STATSFINAL.(textmeasure).nogo_measure1array, 'Sheet1', numtempf);
    xlswrite(stempp,{'Measure1 channel strength NOGO'},'Sheet1', numtempA)
    numtemp1=numtemp1+1+N;
    
    % measure 2
    numtempf=['B' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp,STATSFINAL.(textmeasure).go_measure2array, 'Sheet1', numtempf);
    xlswrite(stempp,{'Measure 2 couple strength GO'},'Sheet1', numtempA)
    N=size(STATSFINAL.(textmeasure).go_measure2array,1)
    numtemp1=numtemp1+1+N;
    
    numtempf=['B' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp,{'Measure 2 couple strength NOGO'},'Sheet1', numtempA)
    xlswrite(stempp,STATSFINAL.(textmeasure).nogo_measure2array, 'Sheet1', numtempf);
%% measure3
% titles - headings
    numtemp1=numtemp1+2;
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
    numtempf=['B' num2str(numtemp1) ':B' num2str(numtemp1)]; 
%     numtempC=['C' num2str(numtemp1) ':C' num2str(numtemp1)];  
%     numtempD=['D' num2str(numtemp1) ':D' num2str(numtemp1)];  
%     numtempE=['E' num2str(numtemp1) ':E' num2str(numtemp1)];  
    xlswrite(stempp, (titles(4)), 'Sheet1', numtempA);
    xlswrite(stempp, {'areas Go'}, 'Sheet1', numtempf);
%     xlswrite(stempp, {'Str Awake'}, 'Sheet1', numtempC);
%     xlswrite(stempp, {'App Awake'}, 'Sheet1', numtempD);
%     xlswrite(stempp, {'App Awake'}, 'Sheet1', numtempE);

% array
    numtemp1=numtemp1+1;
    numtempf=['B' num2str(numtemp1)];
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];
    xlswrite(stempp,STATSFINAL.(textmeasure).go_measure3array,'Sheet1', numtempf);

% p value
    N=length(STATSFINAL.(textmeasure).go_measure3array);
    numtemp1=numtemp1+N+1;
    numtempf=['B' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp,{'p value go'},'Sheet1', numtempA);
    %xlswrite(stempp,STATSFINAL.(textmeasure).measure3p_go, 'Sheet1', numtempf);
    numtemp1=numtemp1+1;
    numtempf=['B' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
    xlswrite(stempp,{'p value nogo'},'Sheet1', numtempA)
    %xlswrite(stempp,STATSFINAL(11).(textmeasure).measure3p_nogo, 'Sheet1', numtempf);
%

%measure3 nogo
    numtemp1=numtemp1+1;
    numtempf=['B' num2str(numtemp1)];  
    numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)]; 
    xlswrite(stempp,{'areas noGo'},'Sheet1', numtempA)
    xlswrite(stempp,STATSFINAL.(textmeasure).nogo_measure3array, 'Sheet1', numtempf);
    
    numtemp1=numtemp1+1+N;
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