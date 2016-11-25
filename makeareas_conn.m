function rightarea=makeareas_conn(Go1_right, Brain_areas, FR,CR, PR, textmeasure, trigger, name, DOI, side)
%     Brain_areas={'FR', 'PR', 'CR', 'CZ', 'OZ', 'FZ'};% {'FL', 'PL', 'CL', 'CZ', 'OZ', 'FZ'}
%     FR=[1,11,12];
%     PR=[4,5,6];
%     CR=[7,9,10];
%     CZ=8;
%     OZ=3;
%     FZ=2;    
% trigger='go', or 'nogo'
% side='Left' or 'Right'
% numRight=[1,2,20:29]; % 12 el.  Cental electrodes excluded
% numLeft=[4, 6:15,17];
% N2=N(numRight);% N2 = { 'FR2'    'FZ2'    'O2'    'PR4'    'PR2'    'PR6'    'CR2'    'CZ2'    'CR6'    'CR4'    'FR6'    'FR4'
% Brain_areas={'FR', 'PR', 'CR', 'CZ', 'OZ', 'FZ'};% {'FL', 'PL', 'CL', 'CZ', 'OZ', 'FZ'}
% FR=[1,11,12];
% PR=[4,5,6];
% CR=[7,9,10];
% CZ=8;
% OZ=3;
% FZ=2;
    
for m=1:length(Brain_areas)
        switch m 
            case 1, stempR=FR;
            case 2, stempR=PR;
            case 3, stempR=CR;
%             case 4, stempR=CZ;
%             case 5, stempR=OZ;
%             case 6, stempR=FZ;
        end
    for k=1:length(Brain_areas)
        switch k 
        case 1, stempC=FR;
        case 2, stempC=PR;
        case 3, stempC=CR;
%         case 4, stempC=CZ;
%         case 5, stempC=OZ;
%         case 6, stempC=FZ;
        end
    sumytemp=0; for ff=1:length(stempR), for gg=1:length(stempC), temp=Go1_right(stempR(ff), stempC(gg)); sumytemp=temp+sumytemp;end, end, sumytemp=sumytemp/(length(stempC)*length(stempR));
    rightarea(m,k)=sumytemp;
    end
end
  
% kalitera na min midenisoume tin pliroforia sto diagwnio, i kathe
% pliroforia dinei tin ididynami tis kathe perioxis
    
%     figure; imagesc(rightarea); set(gca, 'YTickLabel', Brain_areas); set(gca, 'XTickLabel', Brain_areas);
%       set(gca,'Ytick', 1:length(Brain_areas)); set(gca, 'XTick', 1:length(Brain_areas));title([name '-' textmeasure '-' trigger '-' side])
%       figtemp=['Areas-' name '-' textmeasure '-' trigger '-' side];
%       cd(DOI)
%       cd('Pictures')
%       saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
%       clear figtemp
      
      
      
            
         