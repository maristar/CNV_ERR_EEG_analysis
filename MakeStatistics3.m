% Make statistical analysis 07-12-2011 on grandaverages
clear all 
close all

cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program')
cd GO
files=dir('*');
for kk=1:length(files); 
    filenames{kk,:}=files(kk,:).name;
end
ND=kk; clear kk
nchan=input('Number of channels ')
textmeasuresall={'cor','pcor','DTFdelta','DTFtheta','DTFalpha','DTFbeta','DTFgamma'}
for qq=1:length(textmeasuresall)
    textmeasure=textmeasuresall{qq}
    grandaverage=zeros(nchan);
    grandaverage_nogo=zeros(nchan);
    for kkm=3:ND
        cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program')
        cd GO
        %disp(filenames(kkm));
        disp(kkm)
        stemp=([filenames{kkm}]);
        cd(stemp)
        load resultscor
        resultsALL.XYZ=resultscor.XYZ;
        resultsALL.s=resultscor.s;
        switch textmeasure
            case 'pcor', meantemp=squeeze(mean(resultscor.(textmeasure).result2,3));
            case 'cor', meantemp=squeeze(mean(resultscor.(textmeasure).result1,3));
            case {'DTFdelta','DTFtheta','DTFalpha','DTFbeta','DTFgamma'}, meantemp=resultscor.resultsDTF.(textmeasure).meangamma; %for %dtf
        end
        grandaverage=meantemp+grandaverage;
        clear meantemp resultscor
    end
    cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program')
    mkdir(textmeasure)
    cd(textmeasure)
    resultsALL.(textmeasure).go=grandaverage;
    clear kkm
    %% NO GO
    cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program')
    cd NOGO
    files2=dir('*');
    for kk=1:length(files); 
        filenames2{kk,:}=files2(kk,:).name;
    end
    ND=kk; clear kk;
    for kkm=3:ND
        cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program')
        cd NOGO
        disp(kkm)
        stemp=([filenames2{kkm}]);
        cd(stemp)
        load resultscor
        switch textmeasure
            case 'pcor', meantemp=squeeze(mean(resultscor.(textmeasure).result2,3));
            case 'cor', meantemp=squeeze(mean(resultscor.(textmeasure).result1,3));
            case {'DTFdelta','DTFtheta','DTFalpha','DTFbeta','DTFgamma'}, meantemp=resultscor.resultsDTF.(textmeasure).meangamma; %for %dtf
        end
        grandaverage_nogo=meantemp+grandaverage_nogo;
        clear meantemp 
    cd ..   
    end
    cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program')
    cd(textmeasure)
    resultsALL.(textmeasure).nogo=grandaverage_nogo;
    
    %measure1_nogo=squeeze(measure1_nogo);
    XYZ=resultscor.XYZ;   
    N={'FR2', 'CZ1', 'FCZ' 'FL1', 'FL3', 'FL5','CL3', 'CL1', 'CL5', 'PL5', 'PL1', 'PL3','PZC', 'O1', 'PZP', 'OZ', 'O2', 'PR4', 'PR2', 'PR6', 'CR2', 'CZ2', 'CR6', 'CR4', 'FR6', 'FR4'};
    clear resultscor

    cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program')
        
    s=N;
    h=0;
    clear index_pool kj kl 
    index_pool={1:nchan*nchan}; %% to change this to the value of possible combinations
    for kl=1:nchan
        for kj=1:nchan
            if (kl~=kj)
                h=h+1;
                index_pool{h,:}=[num2str(kj) num2str(kl)];
                a=strcmp(index_pool, [num2str(kl),num2str(kj)]);
                onesz=find(a);
                if isempty(onesz)
                    go=squeeze(grandaverage(kl,kj,:));   
                    nogo=squeeze(grandaverage_nogo(kl,kj,:));
                    [p,h] = signrank(go, nogo);
                %[h,p,ci] =ttest2(awake1,sleep1,0.005,'both', 'unequal');
                    significance(kl,kj,:)=p;
                    h_value(kl,kj,:)=h;
                    clear p, h
                end
                clear onesz a 
            end
        end
    end
    for jj=1:nchan;
        h_value(jj,jj)=0;
        significance(jj,jj)=0;
    end
 
    % a way to plot and see what is significant and what is not. 
    for kl=1:nchan
      for kj=1:nchan
         if significance(kl,kj)>0.05
             significance_plot(kl,kj)=1;
            else
             significance_plot(kl,kj)=0;
         end
      end
    end
    
    resultsALL.(textmeasure).signrank_p=significance_plot;
    resultsALL.(textmeasure).signrank_h=h_value;
    cd(textmeasure)
    save resultsALL resultsALL -v7.3
   
    figure; imagesc(h_value);
    set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
    set(gca, 'YTickLabel', s); set(gca, 'XTickLabel', s);
    axis xy; axis tight; colorbar('location','EastOutside')
    
    figure; imagesc(significance_plot);
    set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
    set(gca, 'YTickLabel', s); set(gca, 'XTickLabel', s);
    axis xy; axis tight; colorbar('location','EastOutside')
    title('significance p value')
    
    
end
cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program')
save resultsALL resultsALL -v7.3

%% Now we have the grandaverage and so we go on for the intra-hemispheric statistics
% a) right hemisphere number of channel; 1, 6, 20, 21, 22, 23, 24, 25, 26,
% 27, 28, 29

textmeasure='cor'
numRight=[1, 17:26]; % 12 electrodes
numLeft=[3,4,5, 6, 7, 8,9, 10, 11,12, 14];

load resultsALL

nchan=26
s=resultsALL.s;
XYZ=resultsALL.XYZ;
corgof=zeros(nchan);
corgofL=zeros(nchan);
mcor=resultsALL.(textmeasure).go;

corgof(numRight, numRight)=mcor(numRight,numRight);
corgofL(numLeft, numLeft)=mcor(numLeft,numLeft); 

figure;imagesc(corgof);axis xy; set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
set(gca, 'YTickLabel', s); set(gca, 'XTickLabel', s);axis xy; axis tight; colorbar('location','EastOutside')
figure; imagesc(mcor);axis xy; set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
set(gca, 'YTickLabel', s); set(gca, 'XTickLabel', s);axis xy; axis tight; colorbar('location','EastOutside')
% head plot
plot2dhead_frontal2(corgof, XYZ, s); 


figure;imagesc(corgofL);axis xy; set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
set(gca, 'YTickLabel', s); set(gca, 'XTickLabel', s);axis xy; axis tight; colorbar('location','EastOutside')

% head plot
plot2dhead_frontal2(corgofL, XYZ, s); 

