function rightarea=makeareas_connAP(Go1, Brain_areas, numAnt,numPos, textmeasure, trigger)
%     Brain_areas={'Ant', 'Pos'};% {'FL', 'PL', 'CL', 'CZ', 'OZ', 'FZ'}
%     numAnt=[1,2,3,5:9, 28,29];
%     numPos=[13:23];
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
            case 1, stempA=numAnt;
            case 2, stempP=numPos;
    end
    for k=1:length(Brain_areas)
        switch k 
        case 1, stempA=numAnt;
        case 2, stempP=numPos;
    end
    sumytemp=0; for ff=1:length(numAnt), for gg=1:length(numPos), temp=Go1(numAnt(ff), numPos(gg)); sumytemp=temp+sumytemp;end, end, sumytemp=sumytemp/(length(numAnt)*length(numPos));
    rightarea(m,k)=sumytemp;
    end
end
    
    figure; imagesc(rightarea); set(gca, 'YTickLabel', Brain_areas); set(gca, 'XTickLabel', Brain_areas);
      set(gca,'Ytick', 1:length(Brain_areas)); set(gca, 'XTick', 1:length(Brain_areas));title([textmeasure '-' trigger])
      figtemp=[textmeasure '-' trigger '- Anterior posterior'];
      saveas(gcf, figtemp, 'fig'); saveas(gcf, figtemp, 'jpeg')
            
         