%Simulated dataset to check DTF and maps accuracy, this program checks the
%e-Connectome program.
% 14.11.2011 MLS
data_stim=ones(26, floor(length(data)/55),55);
timedim=length(data_stim);
timeVec2=(1:timedim)./Fs;
a=sin(2*pi*8.*timeVec2);
for kk=1:55; data_stim(1,:,kk)=a;end
b=sin(2*pi*8.*(timeVec2+10))+5*cos(2*pi*50.*timeVec2);
for kk=1:55; data_stim(20,:,kk)=b;end
Fs=250;
% function [DTFa]=DTF_maria_frontal(low_freq, high_freq, p, Fs, data,thismoment, XYZ, B, name, stemp)
cd D:\RIKSHOSPITALET\CNV_RIKS\Programzs\coord
load XYZ
load B
thismoment=datestr(now);
for jj=1:length(thismoment); if ( thismoment(jj)==':' || thismoment(jj)==' '); thismoment(jj)='-';end; end

%[DTF-test]=DTF_maria_frontal(4, 9, 20, Fs, data_stim,thismoment, XYZ, B, name, 'stimulated-dataset8Hz');
% because the function would not run we paste it here. 

low_freq=4
high_freq=10;
p=20
name='test23012012'
nchan=size(data, 1);
num_epochs=size(data,3);

% DTF calculation 
%matlabpool(4)
tic
for jk=1:num_epochs
       ts=squeeze(data(:,:,jk))'; % 4000 x nchan
%        for cc=1:nchan
%            ts1(:,cc)=(ts(:,cc)-mean(ts(:,cc)))/std(ts(:,cc));
%        end
    gamma2 = DTF(ts,low_freq,high_freq,p,Fs);  %% 10 x 10 x 4
    %plot(SBC,'DisplayName','SBC','YDataSource','SBC');figure(gcf); title(['chan ' B(ik)])
    gamma2_all(jk,:,:,:)=gamma2;
    gamma2=[];
end
toc/60
clear jk ts gamma2 

for kk=1:nchan,
       gamma2_all(:,kk,kk,:)=0; %% set zero the diagonal elements 
     end
clear kk    
        
% gamma = number of epochs x nchan x nchan x 4 (sink, source, frequency index)
cd('D:\RIKSHOSPITALET\CNV_RIKS\ANALYZED DATASETS\')
stemp_dir=[thismoment 'DTF-order' num2str(p) '-freqs-' num2str(low_freq) 'to' num2str(high_freq) '5x-noise'];
mkdir(stemp_dir)
cd(stemp_dir)

% make averages for all epochs 
meanA=squeeze(mean(gamma2_all,1)); % nchan x nchan x Nfreq
meanFreq=mean(meanA, 3); %% 10 x 10
gammaFreq=mean(gamma2_all,4); % num_epochs x nchan x nchan 
meangamma=mean(gammaFreq,1);
% Matrix DTF plots 
figure;imagesc(meanFreq);axis xy; 
set(gca, 'YTickLabel', B); set(gca, 'XTickLabel', B);
set(gca,'Ytick', 1:nchan); set(gca, 'XTick', 1:nchan);
colorbar('location','EastOutside');
title([name 'DTF: ' num2str(low_freq) '-' num2str(high_freq) ' Hz']); 
figure_temp=[name(1:end-4) '-matrix-DTF ' num2str(low_freq) '-' num2str(high_freq) '-Hz'];
saveas(gcf, figure_temp, 'fig')
clear source
color_lim=caxis;
clear figure_temp
%%\ New Place for Assymetry index
timeVec2=(1:size(data,2))./Fs;

% save all variables into the allready existing DTF substructure
DTFa.gamma2_all=gamma2_all;
DTFa.gamma_Freq=gammaFreq;
DTFa.gamma_time=meanA;
DTFa.meangamma=squeeze(meangamma);
DTFa.p=p;
DTFa.low_freq=low_freq;
DTFa.high_freq=high_freq;

disp('DTF(20,1)')
disp(meanFreq(20,1))
disp('DTF(1, 20)')
disp(meanFreq(1, 20))
save DTFa DTFa