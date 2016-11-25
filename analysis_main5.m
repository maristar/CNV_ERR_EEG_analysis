%%Analysis of CNV datasets. Maria L. S. 14.11.2011
% this was the analysis_main4
% this is based on the program: DCoffset_removal_20_10_2011.m
% this works on filtered & epoched datasets
thismoment=datestr(now);
for jj=1:length(thismoment); if ( thismoment(jj)==':' || thismoment(jj)==' '); thismoment(jj)='-';end; end
resultscor.now=thismoment;
clear jj 
%% (1) Load the set with the raw dataset
[ALLEEG EEG CURRENTSET ALLCOM]=eeglab;
%cd('D:\RIKSHOSPITALET\CNV RIKS\RAW DATASETS\')
cd('F:\MEDISIN AT UIO NORWAY PC\DATA_DISK_MEDISIN\frontal lobe trauma\CNV RIKS\RAW DATASETS\Epoched\Go')
EEG=pop_loadset; %('101_02_11_2011.set', 'F:\MEDISIN AT UIO NORWAY PC\DATA_DISK_MEDISIN\frontal lobe trauma\CNV RIKS\RAW DATASETS')
[ALLEEG EEG CURRENTSET]=eeg_store(ALLEEG, EEG); % this puts it into a new ALLEEG dataset1
Fs=EEG.srate;

str_description=[EEG.setname thismoment];
for jj=1:length(str_description); if (str_description(jj)==':' || str_description(jj)=='.'); str_description(jj)='-';end; end
disp(str_description)
% mkdir
%cd('D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\')
cd('F:\MEDISIN AT UIO NORWAY PC\DATA_DISK_MEDISIN\frontal lobe trauma\CNV RIKS\ANALYZED_DATASETS\')
mkdir(str_description)
cd(str_description) 
eeglab redraw

%% select channels
numChannel=[7 106 129 62 72 70 75 83 23 24 27 3 124 123 36 35 40 110 104 109 92 97 91 52 51 59];
EEG = pop_select(EEG,  'channel' , numChannel, 'time', [1.4 4.7]);
eeglab redraw
nchan=size(EEG.chanlocs,2);
% define name  of electrodes & positions
for kk=1:nchan; B{kk}=EEG.chanlocs(1,kk).labels; end
N={'FR2', 'CZ1', 'FL1', 'FL3', 'FL5','CL3', 'CL1', 'CL5', 'PL5', 'PL1', 'PL3','PZC', 'O1', 'PZP', 'OZ', 'O2', 'PR4', 'PR2', 'PR6', 'CR2', 'CZ2', 'CR6', 'CR4', 'FR6', 'FR4', 'CZ'};
clear kk
for g=1:nchan; XYZ(g,1)=EEG.chanlocs(1,g).X; XYZ(g,2)=EEG.chanlocs(1,g).Y; XYZ(g,3)=EEG.chanlocs(1,g).Z; end
% clear g  %% it does not work ... 
% define save name
name=str_description(1:21);
EEG = pop_saveset( EEG, 'filename', str_description, 'check', 'on',  ...
    'filepath',['F:\MEDISIN AT UIO NORWAY PC\DATA_DISK_MEDISIN\frontal lobe trauma\CNV RIKS\ANALYZED_DATASETS\' str_description]); % no pop-up D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\'
eeglab redraw
disp('Write down the noisy channels and epochs!!!!')
eegplot(EEG.data) %if we do not wish to scroll we write eegplot('noui', EEG.data)
%% get data
data=double(EEG.data); % 8 x 1425 x 55 meaning nchan x timeduration x num_epochs

p=model_order_maria(data, 1, 20); % function model order -- to see what order is good. 

num_epochs=size(data,3);
result1=zeros(nchan, nchan, num_epochs);
result2=zeros(nchan, nchan, num_epochs);
for k=1:num_epochs
    tempiii=data(1:nchan,:,k)'; %% it wants first the data points, then the number of channels
   result1(:,:,k)=corrcoef(tempiii);
   result2(:,:,k)= partialcorr(tempiii);% partialcorr corrcoef  corrcoef
    a=squeeze(result1(:,:,k));
%       2-d PLOTS     
%     figure; imagesc(a); set(gca,'Ytick', 1:8); set(gca, 'XTick', 1:8); set(gca, 'YTickLabel', s); set(gca, 'XTickLabel', s);
%     axis xy; axis tight; colorbar('location','EastOutside')
%    clear tempiii a 
 end
clear k  tempii

% set diagonal elements to zero
for k=1:num_epochs;
    for jj=1:nchan;
        result1(jj,jj)=0;
        result2(jj,jj)=0;
    end
end
clear jj k a 

textmeasure1='Correlation'; %% NEW
textmeasure2='Partial Correlation';
direct_temp1=[thismoment '-' textmeasure1];
direct_temp2=[thismoment '-' textmeasure2];
mkdir(direct_temp1);
mkdir(direct_temp2);
%% partial correlation
cd(direct_temp2)
pcor_average=mean(result2(:,:,1:end),3);
% matrix plot


figure;imagesc(pcor_average);
set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
set(gca, 'YTickLabel', B); set(gca, 'XTickLabel', B);
axis xy; axis tight; colorbar('location','EastOutside')
title([name(1:end-4) '-' textmeasure2 ': ' 'average']);
figure_temp=[textmeasure2 '- ' 'average' name(1:end-4)]; saveas(gcf, figure_temp, 'fig')

% lines plot
plot2dhead_frontal(pcor_average, XYZ, B); title([name(1:end-4) '-' textmeasure2 ': average' ]);
figure_temp=['Lines-' textmeasure2 '- ' 'average' name(1:end-4)]; 
saveas(gcf, figure_temp, 'fig')

%list of most strong channels
[PCOR]=majorlist(pcor_average, textmeasure2, now, name, nchan, B)
[chan_strength_norm, chan_appearance_norm]=strength_appearance(PCOR.couple_conn, PCOR.couple_conn_values,B);
% input start
% Send the results to excel file 
stempp=['MostCouples-' textmeasure2 '-stats']; % do not delete this!!!!
xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
titles={name};
xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
xlswrite(stempp,B, 'Sheet1', 'B2');
xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');

textmeasure='pcor'
% save to results.. go to RAW\ 
resultscor.(textmeasure).chan_strength_norm_sleep=chan_strength_norm;
resultscor.(textmeasure).chan_appearance_sleep=chan_appearance_norm;
clear chan_appearance_norm chan_strength_norm
clear board board_values                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

cd ..
%% correlation
cd(direct_temp1)
cor_average=mean(result1(:,:,1:end),3);
% matrix plot
figure;imagesc(cor_average);
set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
set(gca, 'YTickLabel', B); set(gca, 'XTickLabel', B);
axis xy; axis tight; colorbar('location','EastOutside')
title([name(1:end-4) '-' textmeasure1 ': ' 'average']);
figure_temp=[textmeasure1 '- ' 'average' name(1:end-4)]; saveas(gcf, figure_temp, 'fig')

% lines plot
plot2dhead_frontal(cor_average, XYZ, B); title([name '-' textmeasure1 ': ''average']);
figure_temp=['Lines-' textmeasure1 '- ' 'average' name(1:end-4)]; 
saveas(gcf, figure_temp, 'fig')
clear figure_temp
[COR]=majorlist(cor_average, textmeasure1, now, name, nchan, B)
[chan_strength_norm, chan_appearance_norm]=strength_appearance(COR.couple_conn, COR.couple_conn_values,B);
% inpout start
% Send the results to excel file 
stempp=['MostCouples-' textmeasure1 '-stats']; % do not delete this!!!!
xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
titles={name};
xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
xlswrite(stempp,B, 'Sheet1', 'B2');
xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');

textmeasure='cor'
% save to results.. go to RAW\ 
resultscor.(textmeasure).chan_strength_norm_sleep=chan_strength_norm;
resultscor.(textmeasure).chan_appearance_sleep=chan_appearance_norm;
clear chan_appearance_norm chan_strength_norm
clear board board_values                         

clear textmeasure1 textmeasure2 tempiii
cd ..
%% DTF
p=10;
cd(['F:\MEDISIN AT UIO NORWAY PC\DATA_DISK_MEDISIN\frontal lobe trauma\CNV RIKS\ANALYZED_DATASETS\' str_description])
[DTFtheta]=DTF_maria_frontal(4, 8, p, Fs, data,thismoment, XYZ, B, name, str_description);
textmeasure='DTFtheta'
[mc_DTFtheta]=majorlist(DTFtheta.meangamma, textmeasure, now, name, nchan, B)
DTFtheta.couple_conn=mc_DTFtheta.couple_conn;
DTFtheta.couple_conn_values=mc_DTFtheta.couple_conn_values;
% appearance-strength
[chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFtheta.couple_conn, DTFtheta.couple_conn_values,B);
stempp=['MostCouples-' textmeasure1 '-stats']; % do not delete this!!!!
xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
titles={name};
xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
xlswrite(stempp,B, 'Sheet1', 'B2');
xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
% save to results.. go to RAW\ 
resultscor.(textmeasure).chan_strength_norm_sleep=chan_strength_norm;
resultscor.(textmeasure).chan_appearance_sleep=chan_appearance_norm;
clear chan_appearance_norm chan_strength_norm
clear board board_values                 
clear mc_DTFtheta
cd ..

textmeasure='DTFdelta';
[DTFdelta]=DTF_maria_frontal(0.1, 4, p, Fs, data,thismoment, XYZ, B, name, str_description);
[mc_DTFdelta]=majorlist(DTFdelta.meangamma, textmeasure, now, name, nchan, B)
DTFdelta.couple_conn=mc_DTFdelta.couple_conn;
DTFdelta.couple_conn_values=mc_DTFdelta.couple_conn_values;
[chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFdelta.couple_conn, DTFdelta.couple_conn_values,B);
resultscor.(textmeasure).chan_strength_norm_sleep=chan_strength_norm;
resultscor.(textmeasure).chan_appearance_sleep=chan_appearance_norm;
% write to excel
stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
titles={name};
xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
xlswrite(stempp,B, 'Sheet1', 'B2');
xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
% save to results.. go to RAW\ 

clear chan_appearance_norm chan_strength_norm
clear board board_values                 
clear mc_DTFdelta
cd ..

textmeasure='DTFalpha';
[DTFalpha]=DTF_maria_frontal(9, 12, p, Fs, data,thismoment, XYZ, B, name, str_description);
[mc_DTFalpha]=majorlist(DTFalpha.meangamma, textmeasure, now, name, nchan, B)
DTFalpha.couple_conn=mc_DTFalpha.couple_conn;
DTFalpha.couple_conn_values=mc_DTFalpha.couple_conn_values;
[chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFalpha.couple_conn, DTFalpha.couple_conn_values,B);
resultscor.(textmeasure).chan_strength_norm_sleep=chan_strength_norm;
resultscor.(textmeasure).chan_appearance_sleep=chan_appearance_norm;
% write to excel
stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
titles={name};
xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
xlswrite(stempp,B, 'Sheet1', 'B2');
xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
clear chan_appearance_norm chan_strength_norm
clear board board_values      
clear mc_DTFalpha
cd ..

textmeasure='DTFbeta';
[DTFbeta]=DTF_maria_frontal(13, 25, p, Fs, data,thismoment, XYZ, B, name, str_description);
[mc_DTFbeta]=majorlist(DTFbeta.meangamma, textmeasure, now, name, nchan, B)
DTFbeta.couple_conn=mc_DTFbeta.couple_conn;
DTFbeta.couple_conn_values=mc_DTFbeta.couple_conn_values;
[chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFbeta.couple_conn, DTFbeta.couple_conn_values,B);
stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
titles={name};
xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
xlswrite(stempp,B, 'Sheet1', 'B2');
xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
% save to results.. go to RAW\ 
resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
resultscor.(textmeasure).chan_appearance=chan_appearance_norm;
clear chan_appearance_norm chan_strength_norm
clear board board_values      
clear mc_DTFbeta
cd ..

textmeasure='DTFgamma';
[DTFgamma]=DTF_maria_frontal(25, 60, p, Fs, data,thismoment, XYZ, B, name, str_description);
[mc_DTFgamma]=majorlist(DTFgamma.meangamma, textmeasure, now, name, nchan, B)
DTFgamma.couple_conn=mc_DTFgamma.couple_conn;
DTFgamma.couple_conn_values=mc_DTFgamma.couple_conn_values;
[chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFgamma.couple_conn, DTFgamma.couple_conn_values,B);
resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
resultscor.(textmeasure).chan_appearance=chan_appearance_norm;

stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
titles={name};
xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
xlswrite(stempp,B, 'Sheet1', 'B2');
xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
% save to results.. go to RAW\ 

clear chan_appearance_norm chan_strength_norm
clear board board_values  
clear mc_DTFgamma
cd ..

%% save the results
resultscor.cor=result1;
resultscor.cor_couples=COR;
resultscor.partcor=result2;
resultscor.pcor_couples=PCOR;
resultscor.nchan=nchan;
resultscor.num_epochs=num_epochs;
resultscor.date=date;
resultscor.epocheddata=data;
resultscor.s=B;
resultscor.XYZ=XYZ;
resultscor.DTFtheta=DTFtheta;
resultscor.DTFdelta=DTFdelta;
resultscor.DTFalpha=DTFalpha;
resultscor.DTFbeta=DTFbeta;
resultscor.DTFgamma=DTFgamma;

save resultscor resultscor -v7.3

% %% wavelet analysis
% width = input('With starting width     ');
% freqN = input('frequency to start?        ');
% repeats = input('For how many times -subsequent frequency bands?   ');
% tic
% %% Define frequencies
% freq1=freqN;
% freqN=freq1+20;  
% step=0.2;
% freqVec =freq1:step:freqN; % 2:0.05:16
% disp(freq1)
% disp(freqN)
% %%
% epoch_length=size(data,2);
% timeVec=(1:length(epoch_length))./Fs;
% Bwav = zeros(length(freqVec), size(data, 2)); %% freqVec x timeLength
% Bwav(1,:)=0;
% %% start!!!
% for jjk=1:nchan,  %size(data10,2)
%     for kk=1:num_epochs
%         temp=data(jjk,:,kk);
%         mean_temp=(temp-mean(temp))/std(temp);
%         % wavelet here!!!
%         for ff=1:length(freqVec)
%             a=(energyvec(freqVec(ff), mean_temp, Fs, width)); 
%             Bwav(ff,:)=a;
%             clear a
%         end
%         Bw_all(:,:,jjk,kk)=Bwav;
%         clear temp mean_temp event_dp kk
%     end
% end
% % epocheddata10 1501 x 399 x 10  (length epoch x num epochs x nchan)
% clear jjk
% 
% %averaging the potentials 
% 
% 
% % avewraging the wavelets per channel
% for jjk=1:nchan
%     chan_wav=squeeze(Bw_all(:,:,jjk,:));    
%     averagechan_wav=mean(chan_wav, 3);
%     averaged_wav10_2(:,:,jjk)=averagechan_wav*1000;
%     clear chan_wav averagechan_wav 
% end
% 
% % plotting the wavelets
% for kkj=1:nchan
%     figure; 
%     temp=squeeze(averaged_wav10_2(:,:, kkj));
%     h=imagesc(timeVec, freqVec, temp); axis xy; title(['channel ' B{kkj}]);xlabel('time (ms)'); ylabel('frequency (Hz)')
%     zlimx(kkj,:,:)=caxis;
%     temp_image=['Averaged wavelet channel ' B{kkj}]
%     saveas(h, temp_image, 'fig')
%     clear temp
% end
% 
% CLim_min=(min(min(squeeze(zlimx))))
% CLim_max=(max(max((zlimx))))
% clims=[CLim_min 600];
% % plotting the wavelets
% for kkj=1:nchan
%     figure; 
%     temp=squeeze(averaged_wav10_2(:,:, kkj));
%     h=imagesc(timeVec, freqVec(12:end), temp(:,12:end), clims); colorbar;axis xy; title(['channel ' B{kkj}]);xlabel('time (ms)'); ylabel('frequency (Hz)')
%    %caxis([CLim_min CLim_max]); 
%     temp_image=['Averaged wavelet channel ' B{kkj}]
%    saveas(h, temp_image, 'fig')
%     clear temp
% end
% 
% figure;
% for kkj=1:nchan
%     temp=squeeze(averaged_wav10_2(:,:, kkj));
%     subplot(2,6,kkj);imagesc(timeVec, freqVec(12:end), temp(:,12:end), clims); axis xy; 
%     xlabel('time (ms)'); ylabel('frequency (Hz)'); title([B{kkj}]);
%    %caxis([CLim_min CLim_max]); 
%     clear temp
% end
% temp_image=['Averaged wavelet channels'];
% saveas(gcf, temp_image, 'fig')
% % plotting the wavelets end
% close all
% FLresultscor.averagedata=averageddata10_2;
% FLresultscor.epocheddata=epocheddata10_2;
% FLresultscor.wavelet=Bw_all;
% FLresultscor.wavelet=averaged_wav10_2;
% % plotting the averages
% figure;
% timeVec=(-start_epoch:length(averageddata10_2)-start_epoch-1);
% for kkk=1:10
%     plot(timeVec, averageddata10_2(:,kkk)); hold on;
% end
%     
% 
% num_epochs=size(epocheddata10_2,2);
% nchan=size(epocheddata10_2,3);
% dataf=epocheddata10_2;
