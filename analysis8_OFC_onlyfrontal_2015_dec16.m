%%Analysis of CNV datasets. Maria L. S. 14.12.2011, corrected 21-02-2012
% this was the analysis_main5
% this is based on the program: DCoffset_removal_20_10_2011.m
% this works on filtered & not epoched datasets, they are epoched here. Cz
% is excluded. 
% additions : new midline frontal area included, data are normalized prior
% to DTF, frequency ranges according to In, 
tic
clear all 
close all
doi='M:\pc\Dokumenter\MATLAB\CNV_dec16\SetFilesFiltered\Interpolate\ICA\ICA_eyeart_removed\'
%'/Users/mstavrin/Documents/MATLAB/CNV/SetFilesFiltered/Interpolate/ICA/' %
% changed 24.11.2014
Todoi='M:\pc\Dokumenter\MATLAB\CNV_dec16\ANALYZED_DATASETS\'
%doi='D:\RIKSHOSPITALET\CNV_RIKS\RAW DATASETS\After_eye_art_removal'% 'D:\OFC\SETS_correct_locations\';
%Todoi='D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED_DATASETS\'%'D:\OFC\ANALYZED_DATASETS'; %

%% Threshold for top 15 strongest couples. 
% crank=15; % Nov 16

%% Look the folders and make the list of the set files in there.
cd(doi)
files=dir('*.set');
for kk=1:length(files); 
    filenames{kk,:}=files(kk,:).name; 
end
clear kk files
% We have 10 folders but appear 12 because the 2 first are . and ..
%filenames= filenames(3:end); 
ND=length(filenames);
%% Define a time stamp
thismoment=datestr(now); 
for jj=1:length(thismoment); if ( thismoment(jj)==':' || thismoment(jj)==' '); thismoment(jj)='-';end; end
clear jj 

% Start 
tic
for kkm=1:ND
    disp(filenames(kkm))
    filename_set_char=char(['Subj_' filenames{kkm}(1:end-4)]);
    disp(kkm)
    cd(Todoi)
    % Make a folder to store the results with the name and the time
    % stamp
    stemp=([filenames{kkm}(1:end-4) '-' thismoment '_only_fr'])
        
    % Correct if the name has dots or empty spaces
    for jj=1:length(stemp); 
        if ( stemp(jj)=='.' || stemp(jj)==' '); stemp(jj)='_';
        end; 
    end
        
    % Make a directory in the analyzed folder and go in. 
    mkdir(stemp)
    cd(stemp)

    % Initialize the resultscor structure to save the data. 
    resultscor.now=thismoment;
    resultscor.name=filenames{kkm}(1:end-4);
    %clear thismoment
    
    %% (1) Load the set with the raw dataset
    [ALLEEG EEG CURRENTSET ALLCOM]=eeglab;
    %cd('D:\RIKSHOSPITALET\CNV RIKS\RAW DATASETS\')
    %cd('F:\MEDISIN AT UIO NORWAY PC\DATA_DISK_MEDISIN\frontal lobe trauma\CNV RIKS\RAW DATASETS\Epoched\Go')
    
    cd(doi)
    
    loadname=[filenames{kkm}] 
    EEG=pop_loadset(loadname); %, STEMP);
    Fs=EEG.srate;
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG)
    %eeglab redraw

    %% Select channels
    Chans_to_take={'E3', 'E123','E124', 'E23', 'E27', 'E24', 'E110', 'E104', 'E109', 'E36', 'E35', 'E40', 'E92', 'E97', 'E91','E51','E52', 'E59', 'E70', 'E75','E83' };
    % Made in 2016 while working with Ingrid
    % Load the names of those channels, names given by us. 
    cd('M:\pc\Dokumenter\MATLAB\CNV\CNV_analysis-master\Newnames_electrodes');
    chan_orig=load('chan_orig.mat')
    for kk=1:length(Chans_to_take)
        temp_chan=Chans_to_take(kk);
        for gg=1:length(chan_orig.chan_orig) % 129
            if strcmp(chan_orig.chan_orig(gg).labels, temp_chan)==1
                numChannel(kk)=gg;
                new_name{kk,:}=chan_orig.chan_orig(gg).newnames
            end
        end
    end
    clear kk gg
    
    % Rename to N 
    N=Chans_to_take;
    clear Chans_to_take;
    %numChannel=[7 11 12 5 6 106 62 72 70 75 83 23 92 97 91 52 51 59]; %numChannel=sort(numChannel)
    EEG = pop_select(EEG,  'channel' , numChannel);
    eeglab redraw
    nchan=size(EEG.chanlocs,2);
    name=stemp(1:21);
    
    % Define just position of electrodes. 
    % Define name  of electrodes & positions
    %         for kk=1:nchan; B{kk}=EEG.chanlocs(1,kk).labels; end
    %         N ={'FR2', 'FZ2','FCZ','CZ1','FZA','FZ1','FL1','FL3','FL5','CL3','CL1','CL5','CR2','CZ2','CR6','CR4','FR6','FR4'};        %N={'FR2', 'CZ1', 'FCZ' 'FL1', 'FL3', 'FL5','CL3', 'CL1', 'CL5', 'PL5', 'PL1', 'PL3','PZC', 'O1', 'PZP', 'OZ', 'O2', 'PR4', 'PR2', 'PR6', 'CR2', 'CZ2', 'CR6', 'CR4', 'FR6', 'FR4'};
    %         clear kk
    for g=1:nchan; 
        XYZ(g,1)=EEG.chanlocs(1,g).X; 
        XYZ(g,2)=EEG.chanlocs(1,g).Y; 
        XYZ(g,3)=EEG.chanlocs(1,g).Z; 
    end
    clear g 

    %% Epoch 
    % Here we can separate in sessions
    sessions={'Go__', 'NoGo'};
    EEG=pop_epoch(EEG, { 'NoGo' }, [-1 4.75], 'newname', [stemp '_sel_epoched']); % Here we change manually, Go__ eller NoGo
    EEG=pop_rmbase(EEG, [-500 -50]);
    eeglab redraw
    % Now cut the useful interval of 1.4 until 4.75
    EEG = pop_select(EEG, 'time', [3.75 4.75]); % new interval 4.10.16 
    % EEG = pop_select(EEG, 'time', [3 4.75]); for the second component 'notrial', [28 34 49 51 56 57]); for 102
    eeglab redraw
    % Note 30.04.2015:
    % Late component 3.2- 3.7 post stimulus (0 at stimulus)
    % Early component 0.5 - 1.0 (post stimulus, 0 at stimulus)
    % end note 30.04.2015
    cd(Todoi)
    cd(stemp)
    [ALLEEG EEG CURRENTSET]=pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', [stemp '_sel_epoched'], 'overwrite', 'on');
    EEG.setname=['Epoched_' stemp];
    eeglab redraw
 

    EEG = pop_saveset(EEG, 'filename', stemp, 'check', 'on')%,  ...
 
    eeglab redraw
    %disp('Write down the noisy channels and epochs!!!!')
    %eegplot(EEG.data) %if we do not wish to scroll we w
    
    %% Get data out of EEGLAB
    data=double(EEG.data); % 8 x 1425 x 55 meaning nchan x timeduration x num_epochs
    
    %% Nan Detection 
    ttt=isnan(data); fff=find(ttt==1);

    %% Model order
    % Clear EEG ALLEEG CURRENTSET CURRENTSTUDY LASTCOM ALLCOM STUDY fff ttt 
    p=model_order_maria(data, 1, 20); % function model order -- to see what order is good. 
   % close (2:21)
    display(['Model order is : ' num2str(p)])
    num_epochs=size(data,3);
    thismoment=resultscor.now;

    %% DTF analysis
    %% Theta band
    [DTFtheta]=DTF_maria_frontal3_tsa(4, 7, p, Fs, data,thismoment, XYZ, N, name, stemp);
    cd(Todoi)
    cd(stemp)
    %cd(['D:\RIKSHOSPITALET\CNV RIKS\ANALYZED DATASETS\' stemp])
    save DTFthetaOct16 DTFtheta
    disp('Done!!!!')
    textmeasure='DTFtheta'
%     
%     [DTFtheta_list]=majorlist_nocrank(DTFtheta.meangamma, textmeasure, now, name, nchan, N, crank)
%     DTFtheta.couple_conn=DTFtheta_list.couple_conn;
%     DTFtheta.couple_conn_values=DTFtheta_list.couple_conn_values;
% % appearance-strength
%     [chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFtheta.couple_conn, DTFtheta.couple_conn_values,N, crank);
%     % br_areas=brainar(chan_strength_norm, chan_appearance_norm);
%     
%     stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
%     xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
%     titles={name};
%     xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
%     xlswrite(stempp,N, 'Sheet1', 'B2');
%     xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
%     xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
%     xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
%     xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
% save to results.. go to RAW\ 
%     resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
%     resultscor.(textmeasure).chan_appearance=chan_appearance_norm;
%     resultscor.(textmeasure).couples=DTFtheta_list; %% na to kanw comment k se ola. to couples yparxei 
    % resultscor.(textmeasure).br_areas=br_areas;
    clear chan_appearance_norm chan_strength_norm
    clear br_areas DTFtheta_list textmeasure
    cd ..
    %% Delta..........
    textmeasure='DTFdelta';
    % Here 4.10.2016
    cd(Todoi)
    cd(stemp)
    [DTFdelta]=DTF_maria_frontal3_tsa(0.1, 4, p, Fs, data,thismoment, XYZ, N, name, stemp);
    cd(Todoi)
    cd(stemp)
    save DTFdelta DTFdelta
    %     [DTFdelta_list]=majorlist_nocrank(DTFdelta.meangamma, textmeasure, now, name, nchan, N,crank)
%     DTFdelta.couple_conn=DTFdelta_list.couple_conn;
%     DTFdelta.couple_conn_values=DTFdelta_list.couple_conn_values;
%     [chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFdelta.couple_conn, DTFdelta.couple_conn_values,N, crank);
%     %br_areas=brainar(chan_strength_norm, chan_appearance_norm);
%     resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
%     resultscor.(textmeasure).chan_appearance_norm=chan_appearance_norm;%%chan_appearance
%     resultscor.(textmeasure).couples=DTFdelta_list;
%     % resultscor.(textmeasure).br_areas=br_areas;
% % write to excel
%     stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
%     xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
%     titles={name};
%     xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
%     xlswrite(stempp,N, 'Sheet1', 'B2');
%     xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
%     xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
%     xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
%     xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
% % save to results.. go to RAW\ 
%     clear chan_appearance_norm chan_strength_norm
%     clear br_areas DTFdelta_list    
    cd ..
%% Alpha
    textmeasure='DTFalpha';
    [DTFalpha]=DTF_maria_frontal3_tsa(8, 13, p, Fs, data,thismoment, XYZ, N, name, stemp);
    cd(Todoi)
    cd(stemp)
    save DTFalpha DTFalpha
    %[DTFalpha_list]=majorlist_nocrank(DTFalpha.meangamma, textmeasure, now, name, nchan, N, crank)
%     DTFalpha.couple_conn=DTFalpha_list.couple_conn;
%     DTFalpha.couple_conn_values=DTFalpha_list.couple_conn_values;
%     [chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFalpha.couple_conn, DTFalpha.couple_conn_values,N, crank);
%     % br_areas=brainar(chan_strength_norm, chan_appearance_norm);
%     resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
%     resultscor.(textmeasure).chan_appearance=chan_appearance_norm;
%     resultscor.(textmeasure).couples=DTFalpha_list;
%     % resultscor.(textmeasure).br_areas=br_areas;
% % write to excel
%     stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
%     xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
%     titles={name};
%     xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
%     xlswrite(stempp,N, 'Sheet1', 'B2');
%     xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
%     xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
%     xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
%     xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
%     clear chan_appearance_norm chan_strength_norm
%     clear DTFbeta_list textmeasure     
%     clear br_areas 
    cd ..
%% Beta
    textmeasure='DTFbeta';
    [DTFbeta]=DTF_maria_frontal3_tsa(14, 30, p, Fs, data,thismoment, XYZ, N, name, stemp);
    cd(Todoi)
    cd(stemp)
    save DTFbeta DTFbeta
%     [DTFbeta_list]=majorlist_nocrank(DTFbeta.meangamma, textmeasure, now, name, nchan, N, crank)
%     DTFbeta.couple_conn=DTFbeta_list.couple_conn;
%     DTFbeta.couple_conn_values=DTFbeta_list.couple_conn_values;
%     [chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFbeta.couple_conn, DTFbeta.couple_conn_values,N, crank);
%     % br_areas=brainar(chan_strength_norm, chan_appearance_norm);
%     stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
%     xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
%     titles={name};
%     xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
%     xlswrite(stempp,N, 'Sheet1', 'B2');
%     xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
%     xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
% 	xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
%     xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
% % save to results.. go to RAW\ 
%     resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
%     resultscor.(textmeasure).chan_appearance=chan_appearance_norm;
%     resultscor.(textmeasure).couples=DTFbeta_list;
%     % resultscor.(textmeasure).br_areas=br_areas;
%     resultscor.N=N;
%     clear chan_appearance_norm chan_strength_norm
%     clear br_areas      
%     clear DTFbeta_list textmeasure
    cd ..
%% Gamma 1
    textmeasure='DTFgamma1';
    
    [DTFgamma1]=DTF_maria_frontal3_tsa(30, 45, p, Fs, data,thismoment, XYZ, N, name, stemp);
    cd(Todoi)
    cd(stemp)
    save DTFgamma1 DTFgamma1
%     [DTFgamma1_list]=majorlist_nocrank(DTFgamma1.meangamma, textmeasure, now, name, nchan, N, crank)
%     DTFgamma1.couple_conn=DTFgamma1_list.couple_conn;
%     DTFgamma1.couple_conn_values=DTFgamma1_list.couple_conn_values;
%     [chan_strength_norm, chan_appearance_norm]=strength_appearance(DTFgamma1.couple_conn, DTFgamma1.couple_conn_values,N, crank);
%     % br_areas=brainar(chan_strength_norm, chan_appearance_norm);
%     resultscor.(textmeasure).chan_strength_norm=chan_strength_norm;
%     resultscor.(textmeasure).chan_appearance_norm=chan_appearance_norm;
%     resultscor.(textmeasure).couples=DTFgamma1_list;
%     % resultscor.(textmeasure).br_areas=br_areas;
%     stempp=['MostCouples-' textmeasure '-stats']; % do not delete this!!!!
%     xlswrite(stempp, {name}, 'Sheet1', 'A1:A1');
%     titles={name};
%     xlswrite(stempp, (titles), 'Sheet1', 'A2:A2');
%     xlswrite(stempp,N, 'Sheet1', 'B2');
%     xlswrite(stempp, {'strength'}, 'Sheet1', 'A3:A3');
%     xlswrite(stempp,chan_strength_norm, 'Sheet1', 'B3');
%     xlswrite(stempp, {'appearance'}, 'Sheet1', 'A4:A4');
%     xlswrite(stempp,chan_appearance_norm, 'Sheet1', 'B4');
% % save to results.. go to RAW\ 
% 
%     clear chan_appearance_norm chan_strength_norm
%     clear board board_values br_areas 
%     
    cd ..

    %% Save the results in a structure 
    cd(Todoi)
    resultscorIngrid.nchan=nchan;
    resultscorIngrid.date=date;
    resultscorIngrid.s=N;
    resultscorIngrid.XYZ=XYZ;
    resultscorIngrid.(filename_set_char).resultsDTF.DTFtheta=DTFtheta;
    resultscorIngrid.(filename_set_char).resultsDTF.DTFdelta=DTFdelta;
    resultscorIngrid.(filename_set_char).resultsDTF.DTFalpha=DTFalpha;
    resultscorIngrid.(filename_set_char).resultsDTF.DTFbeta=DTFbeta;
    resultscorIngrid.(filename_set_char).resultsDTF.DTFgamma1=DTFgamma1;
    
    save resultscorIngrid resultscorIngrid -v7.3
    
    disp(['Max of theta is :' num2str(max(max(DTFtheta.meangamma)))]);
    disp(['Max of delta is :' num2str(max(max(DTFdelta.meangamma)))]);
    disp(['Max of alpha is :' num2str(max(max(DTFalpha.meangamma)))]);
    disp(['Max of beta is :' num2str(max(max(DTFbeta.meangamma)))]);
    disp(['Max of gamma1 is :' num2str(max(max(DTFgamma1.meangamma)))]);
    
close all
clear EEG ALLEEG CURRENTSET data  result1 result2 textmeasure titles ttt stemp stempp num_epochs cor_average cor_list pcor_list resultscor
clear COR PCOR direct_temp1 direct_temp2 p pcor_average fff DTFalpha_list DTFbeta_list DTFdelta_list DTF_gamma1 DTF_gamma1_list DTFtheta_list ans
end
clear kkm
toc/60

%% Make an average array 
% Make an average array of 21 x 21
Average_conn_alpha=zeros(21, 21);
Average_conn_beta=zeros(21, 21);
Average_conn_gamma=zeros(21, 21);
Average_conn_delta=zeros(21, 21);
Average_conn_theta=zeros(21,21);
for kk=1:ND
    disp(filenames(kk))
    filename_set_char=char(['Subj_' filenames{kk}(1:end-4)]);
    % For alpha
    temp=resultscorIngrid.(filename_set_char).resultsDTF.DTFalpha.meangamma;
    Average_conn_alpha=Average_conn_alpha+temp;
    clear temp
    if kk==ND
        Average_conn_alpha=Average_conn_alpha/ND;
    end
    
    % For delta
    temp=resultscorIngrid.(filename_set_char).resultsDTF.DTFdelta.meangamma;
    Average_conn_delta=Average_conn_delta+temp;
    clear temp
    if kk==ND
        Average_conn_delta=Average_conn_delta/ND;
    end
    
    % For Beta
    temp=resultscorIngrid.(filename_set_char).resultsDTF.DTFbeta.meangamma;
    Average_conn_beta=Average_conn_beta+temp;
    clear temp
    if kk==ND
        Average_conn_beta=Average_conn_beta/ND;
    end
    
    % For gamma
    temp=resultscorIngrid.(filename_set_char).resultsDTF.DTFgamma1.meangamma;
    Average_conn_gamma=Average_conn_gamma+temp;
    clear temp
    if kk==ND
        Average_conn_gamma=Average_conn_gamma/ND;
    end
    
    % For theta
    temp=resultscorIngrid.(filename_set_char).resultsDTF.DTFtheta.meangamma;
    Average_conn_theta=Average_conn_theta+temp;
    clear temp
    clear filename_set_char
    if kk==ND
        Average_conn_theta=Average_conn_theta/ND;
    end
end 
clear kk


% Save the results
cd(Todoi)
save Average_conn_all Average_conn*

%% Make the top 15 pairs
crank=15;
N=new_name';
% For theta
textmeasure='theta';
 [DTFtheta_list]=majorlist_all(Average_conn_theta, textmeasure, now, name, nchan, N, crank)

 % For alpha
 textmeasure='alpha';
  [DTFalpha_list]=majorlist_all(Average_conn_alpha, textmeasure, now, name, nchan, N, crank)

% For delta
 textmeasure='delta';
  [DTFdelta_list]=majorlist_all(Average_conn_delta, textmeasure, now, name, nchan, N, crank)

% For gamma
 textmeasure='gamma';
  [DTFgamma_list]=majorlist_all(Average_conn_gamma, textmeasure, now, name, nchan, N, crank)

% For beta
 textmeasure='beta';
  [DTFbeta_list]=majorlist_all(Average_conn_beta, textmeasure, now, name, nchan, N, crank)

 toc 
 
 
 
 
 
 
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