% Make statistical analysis 07-12-2011 on grandaverages
% takes the resultscor and extracts the measure then calculates the
% grandaveage among many channels. Then it calculates the intra-hemispheric
% connectivity. revised 16.12.2011
clear all 
close all

%% GO grandaverage
cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program')
DOI='D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\All interval -new program';
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
        cd(DOI)
        cd GO
        %disp(filenames(kkm));
        disp(kkm)
        stemp=([filenames{kkm}]);
        cd(stemp)
        load resultscor
        resultsALL.XYZ=resultscor.XYZ;
        %resultsALL.s=resultscor.s;
        switch textmeasure
            case 'pcor', meantemp=squeeze(mean(resultscor.(textmeasure).result2,3));
            case 'cor', meantemp=squeeze(mean(resultscor.(textmeasure).result1,3));
            case {'DTFdelta','DTFtheta','DTFalpha','DTFbeta','DTFgamma'}, meantemp=resultscor.resultsDTF.(textmeasure).meangamma; %for %dtf
        end
        grandaverage=meantemp+grandaverage;
        clear meantemp resultscor
    end
    cd(DOI)
    mkdir(textmeasure)
    cd(textmeasure)
    resultsALL.(textmeasure).go=grandaverage;
    clear kkm 
    %% NO GO
    cd(DOI)
    cd NOGO
    files2=dir('*');
    for kk=1:length(files); 
        filenames2{kk,:}=files2(kk,:).name;
    end
    ND=kk; clear kk;
    for kkm=3:ND %% start files for NOGO
        cd(DOI)
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
        clear meantemp resultscor
    cd ..   
    end % end files for Nogo
    cd(DOI)
    cd(textmeasure)
    resultsALL.(textmeasure).nogo=grandaverage_nogo;
    
    %measure1_nogo=squeeze(measure1_nogo);
    XYZ=resultsALL.XYZ;  
    N={'FR2', 'FZ2', 'FCZ', 'CZ1', 'FZA', 'FZ1', 'FL1', 'FL3', 'FL5','CL3', 'CL1', 'CL5', 'PL5', 'PL1', 'PL3','PZC', 'O1', 'PZP', 'OZ', 'O2', 'PR4', 'PR2', 'PR6', 'CR2', 'CZ2', 'CR6', 'CR4', 'FR6', 'FR4'};
    resultsALL.(textmeasure).s=N;
    clear resultscor

    cd(DOI)
       
%     s=N;
%     h=0;
%     clear index_pool kj kl 
%     index_pool={1:nchan*nchan}; %% to change this to the value of possible combinations
%     for kl=1:nchan
%         for kj=1:nchan
%             if (kl~=kj)
%                 h=h+1;
%                 index_pool{h,:}=[num2str(kj) num2str(kl)];
%                 a=strcmp(index_pool, [num2str(kl),num2str(kj)]);
%                 onesz=find(a);
%                 if isempty(onesz)
%                     go=squeeze(grandaverage(kj,kl,:));   
%                     nogo=squeeze(grandaverage_nogo(kj,kl,:));
%                     %[p,h] = signrank(go, nogo);
%                     [h1,p,ci] =ttest2(go,nogo,0.05,'both','unequal');
%                     significance(kl,kj,:)=p;
%                     h_value(kl,kj,:)=h1;
%                     clear p, h1
%                 end % if
%                 clear onesz a 
%             end % if (kl~=kj)
%         end %kj
%     end %kl
%     for jj=1:nchan;
%         h_value(jj,jj)=0;
%         significance(jj,jj)=0;
%     end
%  
%     % a way to plot and see what is significant and what is not. 
%     for kl=1:nchan
%       for kj=1:nchan
%          if significance(kl,kj)<0.05
%              significance_plot(kl,kj)=1;
%             else
%              significance_plot(kl,kj)=0;
%          end
%       end
%     end
%     resultsALL.(textmeasure).ttest=significance;
%     resultsALL.(textmeasure).ttest_plot=significance_plot;
%     resultsALL.(textmeasure).ttest_h=h_value;
%     cd(textmeasure)
%     save resultsALL resultsALL -v7.3
%        
%     figure; imagesc(h_value);title([textmeasure 'h value']);
%     set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
%     set(gca, 'YTickLabel', s); set(gca, 'XTickLabel', s);
%     axis xy; axis tight; colorbar('location','EastOutside')
%     
%     figure; imagesc(significance_plot);title([textmeasure 'p value']);
%     set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
%     set(gca, 'YTickLabel', s); set(gca, 'XTickLabel', s);
%     axis xy; axis tight; colorbar('location','EastOutside')
%     title('significance p value')
%     
%    % clear grandaverage_nogo grandaverage
end
cd(DOI)
save resultsALL resultsALL
%% Now we have the grandaverage and so we go on for the intra-hemispheric statistics
% a) right hemisphere number of channel; 1, 6, 20, 21, 22, 23, 24, 25, 26,

% 27, 28, 29
%% for GO 
s=N;
XYZ=resultsALL.XYZ;
textmeasure='DTFtheta';
mcor=resultsALL.(textmeasure).go;


numRight=[1,2,20:29]; % Cental electrodes excluded
numLeft=[4, 6:15,17];
% s=resultsALL.s;

nchan=length(XYZ);
corgof=zeros(nchan);
corgofL=zeros(nchan);
corgoall=zeros(nchan);
corgof(numRight, numRight)=mcor(numRight,numRight);
corgofL(numLeft, numLeft)=mcor(numLeft,numLeft); 
corgoall=[corgof+corgofL];
thr=0.5*(max(squeeze(max(corgoall))))% to miso tou megistou

grand_lines(textmeasure,corgof,corgofL, XYZ, s, DOI, thr, numLeft, numRight,trigger);

cd(DOI)
save resultsALL resultsALL

%% for NOGO 
s=N;
XYZ=resultsALL.XYZ;
textmeasure='DTFtheta';
mcor=resultsALL.(textmeasure).nogo;


numRight=[1,2,20:29]; % Cental electrodes excluded
numLeft=[4, 6:15,17];
% s=resultsALL.s;

nchan=length(XYZ);
corgof=zeros(nchan);
corgofL=zeros(nchan);
corgoall=zeros(nchan);
corgof(numRight, numRight)=mcor(numRight,numRight);
corgofL(numLeft, numLeft)=mcor(numLeft,numLeft); 
corgoall=[corgof+corgofL];
thr=0.5*(max(squeeze(max(corgoall))))% to miso tou megistou

grand_lines(textmeasure,corgof,corgofL, XYZ, s, DOI, thr, numLeft, numRight,trigger);

cd(DOI)
save resultsALL resultsALL