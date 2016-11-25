% Make statistical analysis 
cd('D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\')
cd GO
files=dir('*');
for kk=1:length(files); 
    filenames{kk,:}=files(kk,:).name;
end
ND=kk;
textmeasuresall={'cor','pcor','DTFdelta','DTFtheta','DTFalpha','DTFbeta','DTFgamma'}
textmeasure='DTFgamma';
nchan=input('Number of channels')
grandaverage=zeros(nchan);
grandaverage_nogo=zeros(nchan);
for kkm=3:ND
    cd('D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\')
    cd GO
    %disp(filenames(kkm));
    disp(kkm)
    stemp=([filenames{kkm}]);
    if stemp~='.'
        cd(stemp)
        load resultscor
        meantemp=squeeze(mean(resultscor.resultsDTF.(textmeasure).meangamma,3));
        grandaverage=meantemp+grandaverage;
%         chan_strength_norm=resultscor.(textmeasure).chan_strength_norm;
%         chan_appearance_norm=resultscor.(textmeasure).chan_appearance_norm;
%          [br_areas]=brainar(chan_strength_norm, chan_appearance_norm);
%          resultscor.(textmeasure).br_areas=br_areas;
%          allcouples=resultscor.(textmeasure).couples.couple_conn_values;
%          save resultscor resultscor -v7.3
        clear meantemp
    end
%     measure1(kkm, :,:)=resultscor.DTFgamma.chan_strength_norm;
%     measure3_areas(kkm,1,:)=resultscor.(textmeasure).br_areas.FL;
%     measure3_areas(kkm,2,:)=resultscor.(textmeasure).br_areas.FR;
%     measure3_areas(kkm,3,:)=resultscor.(textmeasure).br_areas.CZ;
%     measure3_areas(kkm,4,:)=resultscor.(textmeasure).br_areas.CL;
%     measure3_areas(kkm,5,:)=resultscor.(textmeasure).br_areas.CR;
%     measure3_areas(kkm,6,:)=resultscor.(textmeasure).br_areas.PZ;
%     measure3_areas(kkm,7,:)=resultscor.(textmeasure).br_areas.PL;
%     measure3_areas(kkm,8,:)=resultscor.(textmeasure).br_areas.PR;
%     measure3_areas(kkm,9,:)=resultscor.(textmeasure).br_areas.OZ;
    
end
%measure1=squeeze(measure1);

%% NO GO
cd('D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\')
cd NOGO
files2=dir('*');
for kk=1:length(files); 
    filenames2{kk,:}=files2(kk,:).name;
end
ND=kk;

for kkm=3:ND
    cd('D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\')
    cd NOGO
    %disp(filenames2(kkm));
    disp(kkm)
    stemp=([filenames2{kkm}]);
    if stemp~='.'
        cd(stemp)
        load resultscor
        meantemp=squeeze(mean(resultscor.resultsDTF.(textmeasure).meangamma,3));
        grandaverage_nogo=meantemp+grandaverage_nogo;
        %         chan_strength_norm=resultscor.(textmeasure).chan_strength_norm;
%         chan_appearance_norm=resultscor.(textmeasure).chan_appearance_norm;
%          [br_areas]=brainar(chan_strength_norm, chan_appearance_norm);
%          resultscor.(textmeasure).br_areas=br_areas;
%          allcouples_nogo=resultscor.(textmeasure).couples.couple_conn_values;
%          save resultscor resultscor -v7.3
        clear meantemp
        end
%     measure1_nogo(kkm, :,:)=resultscor.DTFgamma.chan_strength_norm;
%     measure3_areas_nogo(kkm,1,:)=resultscor.(textmeasure).br_areas.FL;
%     measure3_areas_nogo(kkm,2,:)=resultscor.(textmeasure).br_areas.FR;
%     measure3_areas_nogo(kkm,3,:)=resultscor.(textmeasure).br_areas.CZ;
%     measure3_areas_nogo(kkm,4,:)=resultscor.(textmeasure).br_areas.CL;
%     measure3_areas_nogo(kkm,5,:)=resultscor.(textmeasure).br_areas.CR;
%     measure3_areas_nogo(kkm,6,:)=resultscor.(textmeasure).br_areas.PZ;
%     measure3_areas_nogo(kkm,7,:)=resultscor.(textmeasure).br_areas.PL;
%     measure3_areas_nogo(kkm,8,:)=resultscor.(textmeasure).br_areas.PR;
%     measure3_areas_nogo(kkm,9,:)=resultscor.(textmeasure).br_areas.OZ;
    cd ..
end
% measure1_nogo=squeeze(measure1_nogo);
XYZ=resultscor.XYZ;   
N={'FR2', 'CZ1', 'FCZ' 'FL1', 'FL3', 'FL5','CL3', 'CL1', 'CL5', 'PL5', 'PL1', 'PL3','PZC', 'O1', 'PZP', 'OZ', 'O2', 'PR4', 'PR2', 'PR6', 'CR2', 'CZ2', 'CR6', 'CR4', 'FR6', 'FR4'};
cd('D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\')
% plot2dhead_frontal(grandaverage, XYZ, N); title(['Grandaverage Go-' textmeasure]);
%     figure_temp=['Lines-' textmeasure '- ' 'Grand_average']; 
%     saveas(gcf, figure_temp, 'fig')
%     
% plot2dhead_frontal(grandaverage_nogo, XYZ, N); title(['Grandaverage NoGo-' textmeasure]);
%     figure_temp=['Lines-' textmeasure '- ' 'Grand_averageNoGo']; 
%     saveas(gcf, figure_temp, 'fig')

figure; plot2deeg3_frontal(grandaverage, XYZ, N); title(['Grandaverage Go-' textmeasure]);
figure_temp=['Grandaverage Go-' textmeasure];
saveas(gcf, figure_temp, 'fig');
    
figure; plot2deeg3_frontal(grandaverage_nogo, XYZ,N); title(['Grandaverage NoGo-' textmeasure]);
figure_temp=['Grandaverage NoGo-' textmeasure];
saveas(gcf, figure_temp, 'fig');

for jj=3:18, [p(jj), h(jj)]=ttest2(measure1(jj,:), measure1_nogo(jj,:)); end  %% file-wise statistics
figure; plot(h, '*');
for kk=1:26, [pc(kk), hc(kk)]=ttest2(measure1(3:end,kk), measure1_nogo(3:end,kk)); end %% channel-wise statistics
figure; plot(hc, '*'); 