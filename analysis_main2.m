%%Analysis of CNV datasets. Maria L. S. 2.11.2011
% this is based on the program: DCoffset_removal_20_10_2011.m

thismoment=datestr(now);
for jj=1:length(thismoment); if ( thismoment(jj)==':' || thismoment(jj)==' '); thismoment(jj)='-';end; end
resultscor.now=thismoment;
clear jj 
%% (1) Load the set with the raw dataset
[ALLEEG EEG CURRENTSET ALLCOM]=eeglab;
cd('F:\MEDISIN AT UIO NORWAY PC\DATA_DISK_MEDISIN\frontal lobe trauma\CNV RIKS\')
EEG=pop_loadset; %('101_02_11_2011.set', 'F:\MEDISIN AT UIO NORWAY PC\DATA_DISK_MEDISIN\frontal lobe trauma\CNV RIKS\RAW DATASETS')
[ALLEEG EEG CURRENTSET]=eeg_store(ALLEEG, EEG); % this puts it into a new ALLEEG dataset1
Fs=EEG.srate;

str_description='connectivity 2';
% mkdir
cd('F:\MEDISIN AT UIO NORWAY PC\DATA_DISK_MEDISIN\frontal lobe trauma\CNV RIKS\')
mkdir(str_description)
cd(str_description)

%% select channels
numChannel=[23 3 36 104 60 85 70 83 50 101];
EEG = pop_select(EEG,  'channel' , [23 3 36 104 60 85 70 83 50 101]);%'point'
eeglab redraw
nchan=size(EEG.chanlocs,2);
% define name  of electrodes & positions
for kk=1:nchan; B{kk}=EEG.chanlocs(1,kk).labels; end
clear kk
for g=1:nchan; XYZ(g,1)=EEG.chanlocs(1,g).X; XYZ(g,2)=EEG.chanlocs(1,g).Y; XYZ(g,3)=EEG.chanlocs(1,g).Z; end
clear g
% define save name
name=EEG.comments(end-24:end);
%% remove DC from all channels - filter data
% Low pass filter 
h=fdesign.lowpass('Fp,Fst,Ap,Ast',0.001,0.05,1,8,250);
d=design(h,'equiripple');
fvtool(d)
data=EEG.data;
nchan=size(data,1);
for k=1:nchan
temp=double(data(k,:));
%temp_filt=filtfilt(d.Numerator,1,temp);
temp_filt=filter(d,temp);
data_filt(k,:)=data(k,:)-temp_filt;
data_filt(k,:)=(data_filt(k,:)-mean(squeeze(data_filt(k,:))));
clear temp temp_filt
end
% check how it looks
ch12=data_filt;
figure;
nn=size(ch12,1);
for kk=1:nn
    a=squeeze(ch12(kk,:));
    plot(a);hold on;
    clear a
end
hold on;
%figure;
raw_ch12=data;
for kk=1:nn
    a=squeeze(raw_ch12(kk,:));
    plot(a, 'g');hold on;legen1=legend('green=raw and blue=filtered');
    clear a
end
clear nn kk k legen1 raw_ch12 ch12 exc d h ans
saveas(gcf,'filteredvsraw', 'fig')
%% put back the data filtered into the EEG structure
EEG.data=data_filt;
EEG = pop_saveset( EEG, 'filename', str_description, 'check', 'on',  ...
    'filepath',['F:\MEDISIN AT UIO NORWAY PC\DATA_DISK_MEDISIN\frontal lobe trauma\CNV RIKS\' str_description]  ); % no pop-up
eeglab redraw
%% Cut first segment
%EEG = pop_select(EEG,  'point' , [30000 240421]);%'point' THis is for the
%iir filters
%% Epoch 
EEG=pop_epoch(EEG, { 'Go__' }, [-1 4.7], 'newname', [str_description '_sel_epoched']);
[ALLEEG EEG CURRENTSET]=pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', [str_description '_sel_epoched'], 'overwrite', 'on');
EEG.setname=['Epoched_' str_description];
EEG=pop_rmbase(EEG, [-1000 -50]);
eeglab redraw

EEG = pop_saveset( EEG,  'filename', [str_description '_sel_epoched_final'],'filepath', ['F:\MEDISIN AT UIO NORWAY PC\DATA_DISK_MEDISIN\frontal lobe trauma\CNV RIKS\' str_description]); %
eeglab redraw
%% apply fft to filtered & epoched data
data=double(EEG.data); % 8 x 1425 x 55 meaning nchan x timeduration x num_epochs
num_epochs=size(data,3)
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
figure;imagesc(pcor_average);title([name(1:end-4) '-' textmeasure2 ': ' 'average']);
set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
set(gca, 'YTickLabel', B); set(gca, 'XTickLabel', B);
axis xy; axis tight; colorbar('location','EastOutside')
figure_temp=[textmeasure2 '- ' 'average' name(1:end-4)]; saveas(gcf, figure_temp, 'fig')
clear figure_temp
% lines plot
figure;
plot2dhead_frontal(pcor_average, XYZ, B); title([name(1:end-4) '-' textmeasure2 ': ''average']);
figure_temp=['Lines-' textmeasure '- ' 'average' name(1:end-4)]; 
saveas(gcf, figure_temp, 'fig')
clear figure_temp
cd ..
%% correlation
cd(direct_temp1)
cor_average=mean(result1(:,:,1:end),3);
% matrix plot
figure;imagesc(cor_average);title([name(1:end-4) '-' textmeasure1 ': ' 'average']);
set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan); 
set(gca, 'YTickLabel', B); set(gca, 'XTickLabel', B);
axis xy; axis tight; colorbar('location','EastOutside')
figure_temp=[textmeasure1 '- ' 'average' name(1:end-4)]; saveas(gcf, figure_temp, 'fig')
clear figure_temp
% lines plot
figure;
plot2dhead_frontal(cor_average, XYZ, B); title([name '-' textmeasure1 ': ''average']);
figure_temp=['Lines-' textmeasure1 '- ' 'average' name(1:end-4)]; 
saveas(gcf, figure_temp, 'fig')
clear figure_temp
%% DTF

 [DTFtheta]=DTF_maria_frontal(4, 8, 20, Fs, data,thismoment, XYZ, B, name, str_description);
 [DTFdelta]=DTF_maria_frontal(0.1, 4, 20, Fs, data,thismoment, XYZ, B, name, str_description);
 [DTFalpha]=DTF_maria_frontal(9, 12, 20, Fs, data,thismoment, XYZ, B, name, str_description);
 [DTFbeta]=DTF_maria_frontal(13, 25, 20, Fs, data,thismoment, XYZ, B, name, str_description);
 [DTFgamma]=DTF_maria_frontal(25, 60, 20, Fs, data,thismoment, XYZ, B, name, str_description);

%% save the results
resultscor.cor=result1;
resultscor.partcor=result2;
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