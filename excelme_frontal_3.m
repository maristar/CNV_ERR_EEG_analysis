function excelme_frontal2(STATSFINAL, Sheet1, textmeasure, thismoment,s)% 
%Send the results to excel file 29-2-2012 Maria Stavrinou
 
ND=size(STATSFINAL(14).(textmeasure).go_measure0array,1);
stempp=['stats' thismoment '-' textmeasure]; 
for kkm=1:(ND)
        name=STATSFINAL(kkm).(textmeasure).name(1:5);
        N=1; 
        numtemp1=0;
        stemppsheet=['Sheet' num2str(kkm)];
        xlswrite(stempp, {thismoment}, stemppsheet, 'A1:A1');
        xlswrite(stempp, {name}, stemppsheet, 'A2:A2');
        
        % Left areas
        % GO
        numtemp1=numtemp1+N+2+1;
        Brain_areas={'FL','CL','PL'};
        numtempC=['C' num2str(numtemp1)];  
        numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
        xlswrite(stempp,{'Go Left 3 areas'},stemppsheet, numtempA)
        xlswrite(stempp, (Brain_areas), stemppsheet, numtempC);

        %arrays
        numtemp1=numtemp1+1;
        numtempB=['B' num2str(numtemp1)]; 
        numtempf=['C' num2str(numtemp1)];
        xlswrite(stempp,(Brain_areas)',stemppsheet, numtempB)
        xlswrite(stempp, STATSFINAL(kkm).(textmeasure).Go_left_areas, stemppsheet, numtempf);
   
        %NoGO
         N=size(STATSFINAL(kkm).(textmeasure).Go_left_areas,1); % 3
         numtemp1=numtemp1+N+2+1;
        %titles
        numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
        numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
        xlswrite(stempp,{'NoGo Left 3 areas'},stemppsheet, numtempA)
        xlswrite(stempp, (Brain_areas), stemppsheet, numtempC);
    
        %arrays
        numtemp1=numtemp1+1;
        numtempB=['B' num2str(numtemp1)]; 
        numtempf=['C' num2str(numtemp1)];
        xlswrite(stempp,(Brain_areas)',stemppsheet, numtempB)
        xlswrite(stempp,STATSFINAL(kkm).(textmeasure).NoGo_left_areas, stemppsheet, numtempf);
        % here 25-02-2012
    
        % RIGHT HEMISPHERe
        % GO
        numtemp1=numtemp1+N+2+1;
        Brain_areas={'FR','CR','PR'};
        numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
        numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
        xlswrite(stempp,{'Go Right 3 areas'},stemppsheet, numtempA)
        xlswrite(stempp, (Brain_areas), stemppsheet, numtempC);

        %arrays
        numtemp1=numtemp1+1;
        numtempB=['B' num2str(numtemp1)]; 
        numtempf=['C' num2str(numtemp1)];
        xlswrite(stempp,(Brain_areas)',stemppsheet, numtempB)
        xlswrite(stempp,STATSFINAL(kkm).(textmeasure).Go_right_areas, stemppsheet, numtempf);

        %NoGO
         N=size(STATSFINAL(kkm).(textmeasure).NoGo_right_areas,1); % 3
         numtemp1=numtemp1+N+2+1;
        %titles
        numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
        numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
        xlswrite(stempp,{'NoGo Right 3 areas'},stemppsheet, numtempA)
        xlswrite(stempp, (Brain_areas),stemppsheet, numtempC);

        %arrays
        numtemp1=numtemp1+1;
        numtempB=['B' num2str(numtemp1)]; 
        numtempf=['C' num2str(numtemp1)];
        xlswrite(stempp,(Brain_areas)',stemppsheet, numtempB)
        xlswrite(stempp,STATSFINAL(kkm).(textmeasure).NoGo_right_areas, stemppsheet, numtempf);
end

name='Grandaverage';
kkm=ND+1;
N=1;
        numtemp1=0;
        stemppsheet=['Sheet' num2str(kkm)];
        xlswrite(stempp, {thismoment}, stemppsheet, 'A1:A1');
        xlswrite(stempp, {name}, stemppsheet, 'A2:A2');
        
        % Left areas
        % GO
        numtemp1=numtemp1+N+2+1;
        Brain_areas={'FL','CL','PL'};
        numtempC=['C' num2str(numtemp1)];  
        numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
        xlswrite(stempp,{'Go Left 3 areas'},stemppsheet, numtempA)
        xlswrite(stempp, (Brain_areas), stemppsheet, numtempC);

        %arrays
        numtemp1=numtemp1+1;
        numtempB=['B' num2str(numtemp1)]; 
        numtempf=['C' num2str(numtemp1)];
        xlswrite(stempp,(Brain_areas)',stemppsheet, numtempB)
        xlswrite(stempp, STATSFINAL(kkm).(textmeasure).Go_left_areas_grandav, stemppsheet, numtempf);
   
        %NoGO
         N=size(STATSFINAL(kkm).(textmeasure).Go_left_areas_grandav,1); % 3
         numtemp1=numtemp1+N+2+1;
        %titles
        numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
        numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
        xlswrite(stempp,{'NoGo Left 3 areas'},stemppsheet, numtempA)
        xlswrite(stempp, (Brain_areas), stemppsheet, numtempC);
    
        %arrays
        numtemp1=numtemp1+1;
        numtempB=['B' num2str(numtemp1)]; 
        numtempf=['C' num2str(numtemp1)];
        xlswrite(stempp,(Brain_areas)',stemppsheet, numtempB)
        xlswrite(stempp,STATSFINAL(kkm).(textmeasure).NoGo_left_areas_grandav, stemppsheet, numtempf);
        % here 25-02-2012
    
        % RIGHT HEMISPHERe
        % GO
        numtemp1=numtemp1+N+2+1;
        Brain_areas={'FR','CR','PR'};
        numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
        numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
        xlswrite(stempp,{'Go Right 3 areas'},stemppsheet, numtempA)
        xlswrite(stempp, (Brain_areas), stemppsheet, numtempC);

        %arrays
        numtemp1=numtemp1+1;
        numtempB=['B' num2str(numtemp1)]; 
        numtempf=['C' num2str(numtemp1)];
        xlswrite(stempp,(Brain_areas)',stemppsheet, numtempB)
        xlswrite(stempp,STATSFINAL(kkm).(textmeasure).Go_right_areas_grandav, stemppsheet, numtempf);

        %NoGO
         N=size(STATSFINAL(kkm).(textmeasure).NoGo_right_areas_grandav,1); % 3
         numtemp1=numtemp1+N+2+1;
        %titles
        numtempC=['C' num2str(numtemp1)];% ':C' num2str(numtemp1)];  
        numtempA=['A' num2str(numtemp1) ':A' num2str(numtemp1)];  
        xlswrite(stempp,{'NoGo Right 3 areas'},stemppsheet, numtempA)
        xlswrite(stempp, (Brain_areas),stemppsheet, numtempC);

        %arrays
        numtemp1=numtemp1+1;
        numtempB=['B' num2str(numtemp1)]; 
        numtempf=['C' num2str(numtemp1)];
        xlswrite(stempp,(Brain_areas)',stemppsheet, numtempB)
        xlswrite(stempp,STATSFINAL(kkm).(textmeasure).NoGo_right_areas_grandav, stemppsheet, numtempf);
       