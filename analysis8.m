%%Analysis of CNV datasets. Maria L. S. 14.12.2011
% this was the analysis_main5
% this is based on the program: DCoffset_removal_20_10_2011.m
% this works on filtered & not epoched datasets, they are epoched here. Cz
% is excluded. 
% additions : new midline frontal area included, data are normalized prior
% to DTF, frequency ranges according to In, 

clear all 
close all
doi='D:\OFC\SETS_correct_locations\'; %'D:\RIKSHOSPITALET\CNV_RIKS\RAW DATASETS\After_eye_art_removal'
Todoi='D:\OFC\ANALYZED_DATASETS'; %'D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\'
cd(doi)
% look the folders and make the list of them
files=dir;
for kk=1:length(files); filenames{kk,:}=files(kk,:).name; end
% We have 10 folders but appear 12 because the 2 first are . and ..
filenames=filenames(3:end); 
ND=length(filenames);
clear kk
thismoment=datestr(now);
for jj=1:length(thismoment); if ( thismoment(jj)==':' || thismoment(jj)==' '); thismoment(jj)='-';end; end
clear jj files
tic
for kkm=1:ND
    disp(filenames(kkm))
    disp(kkm)
        cd(Todoi)
        stemp=([filenames{kkm}(1:end-4) '-' thismoment])
        mkdir(stemp)
        cd(stemp)
        resultscor.now=thismoment;
        resultscor.name=filenames{kkm}(1:end-4);
        %% (1) Load the set with the raw dataset
        [ALLEEG EEG CURRENTSET ALLCOM]=eeglab;
%cd('D:\RIKSHOSPITALET\CNV RIKS\RAW DATASETS\')
%cd('F:\MEDISIN AT UIO NORWAY PC\DATA_DISK_MEDISIN\frontal lobe trauma\CNV RIKS\RAW DATASETS\Epoched\Go')
        STEMP=['D:\OFC\SETS_correct_locations\' filenames{kkm} '_filt']
        EEG=pop_loadset(filenames{kkm}, STEMP);
        Fs=EEG.srate;
        eeglab redraw

        %% select channels
        numChannel=[7 11 12 5 6 106 62 72 70 75 83 23 24 27 3 124 123 36 35 40 110 104 109 92 97 91 52 51 59]; %numChannel=sort(numChannel)
        EEG = pop_select(EEG,  'channel' , numChannel);
        eeglab redraw
        nchan=size(EEG.chanlocs,2);
        name=stemp(1:21);
% define name  of electrodes & positions
        for kk=1:nchan; B{kk}=EEG.chanlocs(1,kk).labels; end
        N ={'FR2', 'FZ2','FCZ','CZ1','FZA','FZ1','FL1','FL3','FL5','CL3','CL1','CL5','PL5','PL1','PL3','PZC','O1','PZP','OZ','O2','PR4','PR2','PR6','CR2','CZ2','CR6','CR4','FR6','FR4'};        %N={'FR2', 'CZ1', 'FCZ' 'FL1', 'FL3', 'FL5','CL3', 'CL1', 'CL5', 'PL5', 'PL1', 'PL3','PZC', 'O1', 'PZP', 'OZ', 'O2', 'PR4', 'PR2', 'PR6', 'CR2', 'CZ2', 'CR6', 'CR4', 'FR6', 'FR4'};
        clear kk
        for g=1:nchan; XYZ(g,1)=EEG.chanlocs(1,g).X; XYZ(g,2)=EEG.chanlocs(1,g).Y; XYZ(g,3)=EEG.chanlocs(1,g).Z; end
        clear g 
%% Filter

        
        %% Epoch 
        EEG=pop_epoch(EEG, { 'NoGo' }, [-1 4.7], 'newname', [stemp '_sel_epoched']);
        EEG=pop_rmbase(EEG, [-1000 -50]);
        eeglab redraw
% now cut the useful interval of 1.4 until 4.75
        EEG = pop_select(EEG, 'time', [3 4.75]); %, 'notrial', [28 34 49 51 56 57]); for 102
        eeglab redraw

        [ALLEEG EEG CURRENTSET]=pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', [stemp '_sel_epoched'], 'overwrite', 'on');
        EEG.setname=['Epoched_' stemp];
        eeglab redraw
        %cd('D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\')
        EEG = pop_saveset(EEG, 'filename', stemp, 'check', 'on',  ...
        'filepath',['D:\OFC\ANALYZED DATASETS\' stemp]); % no pop-up D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\'
        eeglab redraw
        disp('Write down the noisy channels and epochs!!!!')
        eegplot(EEG.data) %if we do not wish to scroll we w
%% get data
        data=double(EEG.data); % 8 x 1425 x 55 meaning nchan x timeduration x num_epochs
%% Nan Detection 
        ttt=isnan(data); fff=find(ttt==1)

        p=model_order_maria(data, 1, 20); % function model order -- to see what order is good. 
        %
        close (2:28);
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
    set(gca, 'YTickLabel', N); set(gca, 'XTickLabel', N);
    axis xy; axis tight; colorbar('location','EastOutside')
    title([name(1:end-4) '-' textmeasure2 ': ' 'average']);
    figure_temp=[textmeasure2 '- ' 'average' name(1:end-4)]; saveas(gcf, figure_temp, 'fig')
% lines plot
    plot2dhead_frontal(pcor_average, XYZ, N); title([name(1:end-4) '-' textmeasure2 ': average' ]);
    figure_temp=['Lines-' textmeasure2 '- ' 'average' name(1:end-4)]; 
    saveas(gcf, figure_temp, 'fig')
%list of most strong couples
    [pcor_list]=majorlist(pcor_average, textmeasure2, now, name, nchan, N)
% list of 26 channels with connectivity strength and appearance
    [chan_strength_norm, chan_appearance_norm]=strength_appearance(pcor_list.couple_conn, pcor_list.couple_conn_values,N);
% input start
br_areas=brainar(chan_strength_norm, chan_appearance_norm);
% Send the results to excel file 
    stempp=['MostCouples-' textmeasure2 '-stats']; % do not delete this!!!!
    xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
    titles={name};
    xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
    xlswrite(stempp,N, 'Sheet1', 'B2');
    xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
    xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
    xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
    xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
    % brain areas
%     xlswrite(stempp, {'brain areas results'}, 'Sheet1', 'A5:A5');
%     xlswrite(stempp,br_areas, 'Sheet1', 'B5');
    %
   
    
    textmeasure='pcor';
% save to results.. \ 
    resultscor.data=data;
    resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
    resultscor.(textmeasure).chan_appearance=chan_appearance_norm;
    resultscor.(textmeasure).result2=result2;
    resultscor.(textmeasure).couples=pcor_list;
    resultscor.(textmeasure).br_areas=br_areas;
   
    clear chan_appearance_norm chan_strength_norm
    clear board board_values br_areas                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    cd ..
%% correlation
    cd(direct_temp1)
    cor_average=mean(result1(:,:,1:end),3);
    % matrix plot
    figure;imagesc(cor_average);
    set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
    set(gca, 'YTickLabel', N); set(gca, 'XTickLabel', N);
    axis xy; axis tight; colorbar('location','EastOutside')
    title([name(1:end-4) '-' textmeasure1 ': ' 'average']);
    figure_temp=[textmeasure1 '- ' 'average' name(1:end-4)]; saveas(gcf, figure_temp, 'fig')
    % lines plot
    plot2dhead_frontal(cor_average, XYZ, N); title([name '-' textmeasure1 ': ''average']);
    figure_temp=['Lines-' textmeasure1 '- ' 'average' name(1:end-4)]; 
    saveas(gcf, figure_temp, 'fig')
    clear figure_temp
    [cor_list]=majorlist(cor_average, textmeasure1, now, name, nchan, N)
    [chan_strength_norm, chan_appearance_norm]=strength_appearance(cor_list.couple_conn, cor_list.couple_conn_values,N);
% inpout start
br_areas=brainar(chan_strength_norm, chan_appearance_norm);
% Send the results to excel file 
    stempp=['MostCouples-' textmeasure1 '-stats']; % do not delete this!!!!
    xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
    titles={name};
    xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
    xlswrite(stempp,N, 'Sheet1', 'B2');
    xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
    xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
    xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
    xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
    textmeasure='cor'
% save to results.. \ 
    resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
    resultscor.(textmeasure).chan_appearance=chan_appearance_norm;
    resultscor.(textmeasure).result1=result1;
    resultscor.(textmeasure).couples=cor_list;
    resultscor.(textmeasure).br_areas=br_areas;
    clear chan_appearance_norm chan_strength_norm
    clear board board_values figure_temp br_areas                  

    clear textmeasure1 textmeasure2 tempiii direct_temp1 direct_temp2 fff
    cd ..
%% DTF
    
    %cd(['D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\' stemp])
    [DTFtheta]=DTF_maria_frontal3_tsa(4, 7, p, Fs, data,thismoment, XYZ, N, name, stemp);
    textmeasure='DTFtheta'
    [DTFtheta_list]=majorlist(DTFtheta.meangamma, textmeasure, now, name, nchan, N)
    DTFtheta.couple_conn=DTFtheta_list.couple_conn;
    DTFtheta.couple_conn_values=DTFtheta_list.couple_conn_values;
% appearance-strength
    [chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFtheta.couple_conn, DTFtheta.couple_conn_values,N);
    br_areas=brainar(chan_strength_norm, chan_appearance_norm);
    
    stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
    xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
    titles={name};
    xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
    xlswrite(stempp,N, 'Sheet1', 'B2');
    xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
    xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
    xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
    xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
% save to results.. go to RAW\ 
    resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
    resultscor.(textmeasure).chan_appearance=chan_appearance_norm;
    resultscor.(textmeasure).couples=DTFtheta_list; %% na to kanw comment k se ola. to couples yparxei 
    resultscor.(textmeasure).br_areas=br_areas;
    clear chan_appearance_norm chan_strength_norm
    clear board board_values br_areas                
    cd ..
% delta..........
    textmeasure='DTFdelta';
    [DTFdelta]=DTF_maria_frontal3_tsa(0.1, 4, p, Fs, data,thismoment, XYZ, N, name, stemp);
    [DTFdelta_list]=majorlist(DTFdelta.meangamma, textmeasure, now, name, nchan, N)
    DTFdelta.couple_conn=DTFdelta_list.couple_conn;
    DTFdelta.couple_conn_values=DTFdelta_list.couple_conn_values;
    [chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFdelta.couple_conn, DTFdelta.couple_conn_values,N);
    br_areas=brainar(chan_strength_norm, chan_appearance_norm);
    resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
    resultscor.(textmeasure).chan_appearance_norm=chan_appearance_norm;%%chan_appearance
    resultscor.(textmeasure).couples=DTFdelta_list;
    resultscor.(textmeasure).br_areas=br_areas;
% write to excel
    stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
    xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
    titles={name};
    xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
    xlswrite(stempp,N, 'Sheet1', 'B2');
    xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
    xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
    xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
    xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
% save to results.. go to RAW\ 
    clear chan_appearance_norm chan_strength_norm
    clear board board_values br_areas                 
    cd ..

    textmeasure='DTFalpha';
    [DTFalpha]=DTF_maria_frontal3_tsa(8, 13, p, Fs, data,thismoment, XYZ, N, name, stemp);
    [DTFalpha_list]=majorlist(DTFalpha.meangamma, textmeasure, now, name, nchan, N)
    DTFalpha.couple_conn=DTFalpha_list.couple_conn;
    DTFalpha.couple_conn_values=DTFalpha_list.couple_conn_values;
    [chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFalpha.couple_conn, DTFalpha.couple_conn_values,N);
    br_areas=brainar(chan_strength_norm, chan_appearance_norm);
    resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
    resultscor.(textmeasure).chan_appearance=chan_appearance_norm;
    resultscor.(textmeasure).couples=DTFalpha_list;
    resultscor.(textmeasure).br_areas=br_areas;
% write to excel
    stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
    xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
    titles={name};
    xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
    xlswrite(stempp,N, 'Sheet1', 'B2');
    xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
    xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
    xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
    xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
    clear chan_appearance_norm chan_strength_norm
    clear board board_values      
    clear mc_DTFalpha br_areas 
    cd ..

    textmeasure='DTFbeta';
    [DTFbeta]=DTF_maria_frontal3_tsa(14, 30, p, Fs, data,thismoment, XYZ, N, name, stemp);
    [DTFbeta_list]=majorlist(DTFbeta.meangamma, textmeasure, now, name, nchan, N)
    DTFbeta.couple_conn=DTFbeta_list.couple_conn;
    DTFbeta.couple_conn_values=DTFbeta_list.couple_conn_values;
    [chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFbeta.couple_conn, DTFbeta.couple_conn_values,N);
    br_areas=brainar(chan_strength_norm, chan_appearance_norm);
    stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
    xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
    titles={name};
    xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
    xlswrite(stempp,N, 'Sheet1', 'B2');
    xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
    xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
	xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
    xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
% save to results.. go to RAW\ 
    resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
    resultscor.(textmeasure).chan_appearance=chan_appearance_norm;
    resultscor.(textmeasure).couples=DTFbeta_list;
    resultscor.(textmeasure).br_areas=br_areas;
    resultscor.N=N;
    clear chan_appearance_norm chan_strength_norm
    clear board board_values br_areas      
    clear mc_DTFbeta
    cd ..

    textmeasure='DTFgamma1';
    [DTFgamma1]=DTF_maria_frontal3_tsa(30, 45, p, Fs, data,thismoment, XYZ, N, name, stemp);
    [DTFgamma1_list]=majorlist(DTFgamma1.meangamma, textmeasure, now, name, nchan, N)
    DTFgamma1.couple_conn=DTFgamma1_list.couple_conn;
    DTFgamma1.couple_conn_values=DTFgamma1_list.couple_conn_values;
    [chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFgamma1.couple_conn, DTFgamma.couple_conn_values,N);
    br_areas=brainar(chan_strength_norm, chan_appearance_norm);
    resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
    resultscor.(textmeasure).chan_appearance_norm=chan_appearance_norm;
    resultscor.(textmeasure).couples=DTFgamma1_list;
    resultscor.(textmeasure).br_areas=br_areas;
    stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
    xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
    titles={name};
    xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
    xlswrite(stempp,N, 'Sheet1', 'B2');
    xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
    xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
    xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
    xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
% save to results.. go to RAW\ 

    clear chan_appearance_norm chan_strength_norm
    clear board board_values br_areas 
    
    cd ..

%% save the results
   
    resultscor.nchan=nchan;
    resultscor.num_epochs=num_epochs;
    resultscor.date=date;
    resultscor.epocheddata=data;
    resultscor.s=N;
    resultscor.XYZ=XYZ;
    resultscor.resultsDTF.DTFtheta=DTFtheta;
    resultscor.resultsDTF.DTFdelta=DTFdelta;
    resultscor.resultsDTF.DTFalpha=DTFalpha;
    resultscor.resultsDTF.DTFbeta=DTFbeta;
    resultscor.resultsDTF.DTFgamma1=DTFgamma1;
    save resultscor resultscor -v7.3

close all
clear EEG ALLEEG CURRENTSET data  result1 result2 textmeasure titles ttt stemp stempp num_epochs cor_average
clear COR PCOR DTFtheta DTFdelta DTFalpha DTFbeta DTFgamma direct_temp1 direct_temp2 p pcor_average fff
end
toc/60
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
